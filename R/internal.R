# ---------------------------------------------------------------------------
# Operator utilities
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x


# ---------------------------------------------------------------------------
# Python module access
# ---------------------------------------------------------------------------

#' Retrieve (and lazily initialise) the cached esfs Python module
#'
#' On first call this imports esfs, suppresses the verbose startup banner,
#' sets the matplotlib backend to non-interactive Agg, and caches everything
#' in `.esfs_pkg_env`. Subsequent calls return the cached module instantly.
#'
#' @keywords internal
.get_esfs <- function(conda_env = NULL) {
  # Return cached module immediately if already initialised
  if (!is.null(.esfs_pkg_env$esfs)) return(invisible(.esfs_pkg_env$esfs))

  # --- Environment configuration ---
  # Must happen before Python starts.  If Python is already running we can
  # only warn — the env was (or should have been) set via .onLoad / .Rprofile.
  if (!is.null(conda_env)) {
    if (reticulate::py_available(initialize = FALSE)) {
      warning(
        "Python is already initialised; the 'conda_env' argument has no effect.\n",
        "To use a specific environment, set options(esfs.conda_env = \"...\") ",
        "in your ~/.Rprofile before loading any package.",
        call. = FALSE
      )
    } else {
      if (reticulate::condaenv_exists(conda_env)) {
        reticulate::use_condaenv(conda_env, required = TRUE)
      } else if (reticulate::virtualenv_exists(conda_env)) {
        reticulate::use_virtualenv(conda_env, required = TRUE)
      } else {
        stop(
          "Environment '", conda_env, "' does not exist.\n",
          "Run install_esfs(envname = \"", conda_env, "\") to create it.",
          call. = FALSE
        )
      }
    }
  }

  # --- Importability check with actionable error ---
  if (!reticulate::py_module_available("esfs")) {
    env <- getOption("esfs.conda_env",
                     Sys.getenv("ESFS_CONDA_ENV", unset = "r-esfs"))
    stop(
      "The 'esfs' Python package is not available in the current environment.\n\n",
      "To fix, run:\n",
      "  install_esfs()   # creates '", env, "' environment with esfs\n\n",
      "If esfs is installed in a different environment, add to ~/.Rprofile:\n",
      "  options(esfs.conda_env = \"your-env-name\")\n",
      "then restart R and reload SeuratESFS.\n\n",
      "Run check_esfs() for a full diagnostic report.",
      call. = FALSE
    )
  }

  # --- Import and configure ---
  esfs <- reticulate::import("esfs", convert = FALSE)
  # Suppress the ASCII banner that configure() prints on first import
  esfs$configure(gpu = TRUE, upcast = FALSE, verbose = FALSE)
  .esfs_pkg_env$esfs <- esfs

  # Non-interactive matplotlib backend so figures render to file, not screen
  plt <- reticulate::import("matplotlib.pyplot", convert = FALSE)
  plt$switch_backend("Agg")
  .esfs_pkg_env$plt <- plt

  invisible(esfs)
}


# ---------------------------------------------------------------------------
# Seurat <-> AnnData bridge
# ---------------------------------------------------------------------------

#' Build an AnnData object from a Seurat object
#'
#' Extracts the count matrix (genes × cells) from Seurat, transposes it to
#' (cells × genes), and constructs a Python AnnData.  Cell and feature
#' metadata are attached as `obs` and `var` data frames.
#'
#' @param object A Seurat object.
#' @param assay  Assay to pull counts from (default: `DefaultAssay(object)`).
#' @param layer  Layer name (default: `"counts"`).
#'
#' @return A Python AnnData object (not converted — stays as Python reference).
#' @keywords internal
.build_adata <- function(object,
                         assay = Seurat::DefaultAssay(object),
                         layer = "counts") {
  anndata <- reticulate::import("anndata", convert = FALSE)
  scipy   <- reticulate::import("scipy",   convert = FALSE)

  # Extract sparse count matrix (genes × cells) and transpose to (cells × genes)
  counts <- Seurat::LayerData(object, assay = assay, layer = layer)
  counts_T <- Matrix::t(counts)   # now cells × genes

  # Convert dgCMatrix -> scipy csr_matrix.
  # reticulate handles the dgCMatrix -> scipy conversion automatically.
  X_scipy <- scipy$sparse$csr_matrix(counts_T)

  # Cell and feature metadata as plain data frames
  obs_df <- as.data.frame(object[[]])                         # cell metadata
  var_df <- as.data.frame(object[[assay]][[]])                # gene metadata

  adata <- anndata$AnnData(
    X   = X_scipy,
    obs = obs_df,
    var = var_df
  )

  adata
}


#' Restore previously computed ESFS matrices into an AnnData object
#'
#' Reads ES matrices (ESSs, EPs, SWs, SGs) and the scaled counts layer from
#' the Seurat `misc` slot and populates the corresponding AnnData slots so
#' that downstream ESFS functions see them exactly as if the earlier steps had
#' been run in the same Python session.
#'
#' @param adata  Python AnnData object (modified in place on the Python side).
#' @param object Seurat object holding previously stored results.
#' @param label  Secondary features label (default `"Self"`).
#'
#' @return The AnnData object (for chaining).
#' @keywords internal
.restore_esfs_to_adata <- function(adata, object, label = "Self") {
  np <- reticulate::import("numpy", convert = FALSE)

  # Scaled counts layer
  scaled <- Seurat::Misc(object, slot = "esfs_Scaled_Counts")
  if (!is.null(scaled)) {
    # stored as genes×cells; AnnData wants cells×genes
    adata$layers["Scaled_Counts"] <- np$array(
      Matrix::t(scaled), dtype = "float32"
    )
  }

  # Gene-level ES matrices stored in varm (genes × genes square)
  for (key in c("ESSs", "EPs", "SWs", "SGs")) {
    slot_name <- paste0("esfs_", key)
    mat       <- Seurat::Misc(object, slot = slot_name)
    if (!is.null(mat)) {
      adata$varm[paste0(label, "_", key)] <- np$array(mat, dtype = "float32")
    }
  }

  # Gene weights / ranks from ES_rank_genes (stored in var)
  weights <- Seurat::Misc(object, slot = "esfs_Gene_Weights")
  ranks   <- Seurat::Misc(object, slot = "esfs_Gene_Ranks")
  if (!is.null(weights)) {
    py_weights <- reticulate::r_to_py(as.data.frame(weights))
    adata$var["ESFS_Gene_Weights"] <- py_weights["ESFS_Gene_Weights"]
  }
  if (!is.null(ranks)) {
    py_ranks <- reticulate::r_to_py(as.data.frame(ranks))
    adata$var["ES_Rank"] <- py_ranks["ES_Rank"]
  }

  # KNN distance graph (sparse cells × cells)
  knn_graph <- object[["esfs.knn"]]
  if (!is.null(knn_graph)) {
    scipy <- reticulate::import("scipy", convert = FALSE)
    knn_r <- as(knn_graph, "CsparseMatrix")
    adata$obsp["correlation_distance_kNN"] <- scipy$sparse$csr_matrix(knn_r)
  }

  adata
}


#' Store ESFS results from an AnnData object back into Seurat misc
#'
#' @param object Seurat object to update.
#' @param adata  Python AnnData with newly computed results.
#' @param label  Secondary features label (default `"Self"`).
#'
#' @return Updated Seurat object.
#' @keywords internal
.store_esfs_from_adata <- function(object, adata, label = "Self") {
  # Scaled counts layer (genes × cells in Seurat convention)
  py_keys <- reticulate::py_to_r(reticulate::py_eval(
    "list(adata.layers.keys())", convert = TRUE,
    local = list(adata = adata)
  ))
  if ("Scaled_Counts" %in% py_keys) {
    scaled_py <- adata$layers["Scaled_Counts"]
    scaled_r  <- reticulate::py_to_r(scaled_py)   # cells × genes
    Seurat::Misc(object, slot = "esfs_Scaled_Counts") <- Matrix::t(
      Matrix::Matrix(scaled_r, sparse = FALSE)
    )
  }

  # Gene-level ES matrices (varm stores genes × n)
  varm_keys <- reticulate::py_to_r(reticulate::py_eval(
    "list(adata.varm.keys())", convert = TRUE,
    local = list(adata = adata)
  ))
  for (key in c("ESSs", "EPs", "SWs", "SGs")) {
    varm_key <- paste0(label, "_", key)
    if (varm_key %in% varm_keys) {
      mat <- reticulate::py_to_r(adata$varm[varm_key])
      Seurat::Misc(object, slot = paste0("esfs_", key)) <- mat
    }
  }

  # Gene weights / ranks written by ES_rank_genes into adata.var
  var_df <- reticulate::py_to_r(adata$var)
  if ("ESFS_Gene_Weights" %in% colnames(var_df)) {
    Seurat::Misc(object, slot = "esfs_Gene_Weights") <-
      var_df["ESFS_Gene_Weights"]
  }
  if ("ES_Rank" %in% colnames(var_df)) {
    Seurat::Misc(object, slot = "esfs_Gene_Ranks") <- var_df["ES_Rank"]
  }

  object
}


#' Build AnnData and restore all previously stored ESFS results
#'
#' Convenience wrapper combining `.build_adata()` and
#' `.restore_esfs_to_adata()`.  Every exported wrapper function calls this to
#' get an AnnData that looks like an uninterrupted Python session.
#'
#' @param object Seurat object.
#' @param assay  Assay name.
#' @param label  Secondary features label.
#'
#' @return Python AnnData with previously computed ESFS results repopulated.
#' @keywords internal
.build_and_restore_adata <- function(object,
                                     assay = Seurat::DefaultAssay(object),
                                     label = "Self") {
  adata <- .build_adata(object, assay = assay)
  .restore_esfs_to_adata(adata, object, label = label)
  adata
}


# ---------------------------------------------------------------------------
# Matplotlib display helper
# ---------------------------------------------------------------------------

#' Save the current matplotlib figure to a temp PNG and display it in R
#'
#' Called immediately after any Python plotting function that creates a
#' matplotlib figure as a side-effect.  The figure is saved to a temporary
#' file, loaded into R, and rendered via `grid::grid.raster()`.
#'
#' @param dpi Resolution of the saved figure (default 150).
#'
#' @return Invisibly returns the path to the temp PNG file.
#' @keywords internal
.display_matplotlib_figure <- function(dpi = 150L) {
  plt <- .esfs_pkg_env$plt
  if (is.null(plt)) {
    plt <- reticulate::import("matplotlib.pyplot", convert = FALSE)
  }

  tmp <- tempfile(fileext = ".png")
  plt$savefig(tmp, bbox_inches = "tight", dpi = as.integer(dpi))
  plt$close("all")

  img <- png::readPNG(tmp)
  grid::grid.newpage()
  grid::grid.raster(img)

  invisible(tmp)
}
