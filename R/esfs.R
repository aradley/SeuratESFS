# ===========================================================================
# Preprocessing
# ===========================================================================

#' Scale and filter the count matrix for ESFS
#'
#' Wraps `esfs.create_scaled_matrix()`. Scales gene expression to \[0, 1\]
#' using per-gene percentile clipping and removes genes with low total
#' expression. The scaled matrix is stored in the Seurat object's `misc` slot
#' and restored automatically into AnnData on subsequent ESFS calls.
#'
#' @param object          A Seurat object with raw counts.
#' @param assay           Assay to use (default: `DefaultAssay(object)`).
#' @param clip_percentile Upper percentile used for scaling (default `97.5`).
#' @param log_scale       If `TRUE`, log-transform before scaling.
#' @param min_total_expr  Genes with total expression below this threshold are
#'   removed (default `50`).
#' @param conda_env       Name of the conda/virtualenv containing esfs.
#'   Usually not needed if `options(esfs.conda_env)` is set.
#'
#' @return The Seurat object with scaled counts stored in `misc`.
#' @export
#' @seealso [CalcESMatrices()]
#' @references Radley A. et al. (2026) \url{https://www.biorxiv.org/content/10.64898/2026.01.26.701684v1}
CreateESFSScaledMatrix <- function(
  object,
  assay           = Seurat::DefaultAssay(object),
  clip_percentile = 97.5,
  log_scale       = FALSE,
  min_total_expr  = 50L,
  conda_env       = NULL
) {
  esfs  <- .get_esfs(conda_env)
  adata <- .build_adata(object, assay = assay)

  adata <- esfs$create_scaled_matrix(
    adata,
    clip_percentile       = clip_percentile,
    log_scale             = log_scale,
    Min_Total_Expression  = as.integer(min_total_expr)
  )

  object <- .store_esfs_from_adata(object, adata)
  object <- Seurat::LogSeuratCommand(object)
  object
}


#' Calculate pairwise ES matrices (ES-GSS)
#'
#' Wraps `esfs.parallel_calc_es_matrices()`. Computes four pairwise
#' information-theoretic metrics (ESS, EP, SW, SG) between all genes using
#' the scaled count matrix. This is the core ES-GSS step.
#'
#' GPU acceleration is used automatically when available:
#' - **NVIDIA CUDA**: full GPU computation via CuPy
#' - **Apple Silicon (MLX/Metal)**: GPU-accelerated ES matrix calculations
#' - **CPU**: Numba JIT-compiled parallel computation
#'
#' @param object          A Seurat object that has been processed with
#'   [CreateESFSScaledMatrix()].
#' @param assay           Assay to use.
#' @param label           Label for the secondary features (default `"Self"`
#'   compares each gene to every other gene).
#' @param save_matrices   Which ES matrices to save in the Seurat object.
#'   Default `c("ESSs", "EPs")`; can also include `"SWs"` and `"SGs"`.
#' @param n_cores         Number of CPU cores for the CPU path (`-1L` = all
#'   available). Ignored when GPU is active.
#' @param conda_env       Conda/virtualenv environment name.
#'
#' @return The Seurat object with ES matrices stored in `misc`.
#' @export
#' @seealso [CreateESFSScaledMatrix()], [ESRankGenes()]
CalcESMatrices <- function(
  object,
  assay         = Seurat::DefaultAssay(object),
  label         = "Self",
  save_matrices = c("ESSs", "EPs"),
  n_cores       = -1L,
  conda_env     = NULL
) {
  esfs  <- .get_esfs(conda_env)
  adata <- .build_and_restore_adata(object, assay = assay, label = label)

  adata <- esfs$parallel_calc_es_matrices(
    adata,
    secondary_features_label = label,
    save_matrices            = reticulate::tuple(save_matrices),
    use_cores                = as.integer(n_cores)
  )

  object <- .store_esfs_from_adata(object, adata, label = label)
  object <- Seurat::LogSeuratCommand(object)
  object
}


# ===========================================================================
# Core algorithms
# ===========================================================================

#' Find optimal cluster combinations per gene (ES-CCF)
#'
#' Wraps `esfs.ES_CCF()`. For each gene, identifies the combination of
#' clusters that maximises ES correlation. Requires an over-clustered
#' one-hot-encoded cluster matrix as secondary features in `adata.obsm`.
#'
#' @param object    A Seurat object with ES matrices computed.
#' @param assay     Assay to use.
#' @param label     Secondary features label matching the clustering stored in
#'   the Seurat object.
#' @param n_cores   Number of CPU cores (`-1L` = all).
#' @param conda_env Conda/virtualenv environment name.
#'
#' @return The Seurat object with CCF results stored in `misc`.
#' @export
#' @seealso [CalcESMatrices()], [RunESFMG()]
RunESCCF <- function(
  object,
  assay     = Seurat::DefaultAssay(object),
  label,
  n_cores   = -1L,
  conda_env = NULL
) {
  esfs  <- .get_esfs(conda_env)
  adata <- .build_and_restore_adata(object, assay = assay, label = label)

  adata <- esfs$ES_CCF(
    adata,
    secondary_features_label = label,
    use_cores                = as.integer(n_cores)
  )

  object <- .store_esfs_from_adata(object, adata, label = label)
  object <- Seurat::LogSeuratCommand(object)
  object
}


#' Select N non-redundant marker genes (ES-FMG)
#'
#' Wraps `esfs.ES_FMG()`. Uses simulated annealing to select a set of N genes
#' that collectively capture maximally distinct expression patterns.
#'
#' @param object      A Seurat object with ES matrices computed.
#' @param N           Number of marker genes to select.
#' @param assay       Assay to use.
#' @param label       Secondary features label (default `"Self"`).
#' @param input_genes Optional character vector restricting the candidate gene
#'   pool. `NULL` uses all genes in the ES matrices.
#' @param num_reheats Number of simulated-annealing reheats (default `3`).
#' @param resolution  Redundancy penalty (default `1`). Lower values allow
#'   more similar genes; use `0.5`–`0.7` for gradient/trajectory data.
#' @param conda_env   Conda/virtualenv environment name.
#'
#' @return A named list with elements:
#'   - `genes`: character vector of selected marker gene names
#'   - `scores`: numeric vector of per-gene ES-FMG scores
#'   - `object`: updated Seurat object with marker list in `misc`
#' @export
#' @seealso [CalcESMatrices()], [RunESCCF()]
RunESFMG <- function(
  object,
  N,
  assay       = Seurat::DefaultAssay(object),
  label       = "Self",
  input_genes = NULL,
  num_reheats = 3L,
  resolution  = 1.0,
  conda_env   = NULL
) {
  esfs  <- .get_esfs(conda_env)
  adata <- .build_and_restore_adata(object, assay = assay, label = label)

  py_input <- if (is.null(input_genes)) {
    reticulate::py_none()
  } else {
    reticulate::r_to_py(as.list(input_genes))
  }

  result <- esfs$ES_FMG(
    adata,
    N                        = as.integer(N),
    secondary_features_label = label,
    input_genes              = py_input,
    num_reheats              = as.integer(num_reheats),
    resolution               = resolution
  )

  genes  <- reticulate::py_to_r(result[[1]])
  scores <- reticulate::py_to_r(result[[2]])

  Seurat::Misc(object, slot = "esfs_markers") <-
    data.frame(gene = genes, score = scores, stringsAsFactors = FALSE)

  object <- Seurat::LogSeuratCommand(object)
  list(genes = genes, scores = scores, object = object)
}


# ===========================================================================
# Analysis & ranking
# ===========================================================================

#' Rank genes by ES network connectivity
#'
#' Wraps `esfs.ES_rank_genes()`. Builds a weighted gene network from the ES
#' matrices and ranks genes by normalised network connectivity. Rankings are
#' stored in the Seurat object's `misc` slot.
#'
#' @param object            A Seurat object with ES matrices computed.
#' @param assay             Assay to use.
#' @param label             Secondary features label (default `"Self"`).
#' @param ep_threshold      Minimum EP value for an edge to be included
#'   (default `0`).
#' @param ess_threshold     Minimum ESS value for an edge to be included
#'   (default `0.01`).
#' @param exclude_genes     Character vector of gene names to exclude.
#' @param known_genes       Character vector of genes whose ranks will be
#'   printed to the console (informational only).
#' @param min_edges         Genes with fewer than this many edges are pruned
#'   from the network (default `5`).
#' @param conda_env         Conda/virtualenv environment name.
#'
#' @return The Seurat object with gene weights and ranks stored in `misc`.
#' @export
#' @seealso [CalcESMatrices()], [PlotESFSGeneUMAP()]
ESRankGenes <- function(
  object,
  assay         = Seurat::DefaultAssay(object),
  label         = "Self",
  ep_threshold  = 0.0,
  ess_threshold = 0.01,
  exclude_genes = NULL,
  known_genes   = NULL,
  min_edges     = 5L,
  conda_env     = NULL
) {
  esfs  <- .get_esfs(conda_env)
  adata <- .build_and_restore_adata(object, assay = assay, label = label)

  np <- reticulate::import("numpy", convert = FALSE)

  py_exclude <- if (is.null(exclude_genes)) {
    reticulate::py_none()
  } else {
    reticulate::tuple(as.list(exclude_genes))
  }
  py_known <- if (is.null(known_genes)) {
    np$array(character(0))
  } else {
    np$array(known_genes)
  }

  adata <- esfs$ES_rank_genes(
    adata,
    EP_threshold             = ep_threshold,
    ESS_threshold            = ess_threshold,
    exclude_genes            = py_exclude,
    known_important_genes    = py_known,
    secondary_features_label = label,
    min_edges                = as.integer(min_edges)
  )

  object <- .store_esfs_from_adata(object, adata, label = label)
  object <- Seurat::LogSeuratCommand(object)
  object
}


#' Smooth gene expression by KNN averaging
#'
#' Wraps `esfs.knn_smooth_gene_expression()`. Smooths expression by averaging
#' each cell over its k nearest neighbours in the correlation distance space
#' defined by `use_genes`.
#'
#' **Backend notes:**
#' - KNN computation is GPU-accelerated on NVIDIA (cuML). Apple Silicon (MLX)
#'   and CPU use chunked BLAS matrix multiplication.
#' - The KNN distance graph is always saved to the Seurat Graphs slot
#'   (`"esfs.knn"`) and restored on subsequent calls to avoid recomputation.
#' - The smoothed expression matrix is dense (cells × genes). For large
#'   datasets it may exceed `max_size_gb` and will not be transferred to R.
#'   In that case use the AnnData directly in Python.
#'
#' @param object           A Seurat object.
#' @param use_genes        Character vector of gene names used to compute
#'   cell-cell distances.
#' @param assay            Assay to use.
#' @param knn              Number of nearest neighbours (default `30`).
#' @param log_scale        Log2-transform expression before distance
#'   computation (default `FALSE`).
#' @param force_recalculate If `TRUE`, recompute KNN even if a stored graph
#'   exists (default `FALSE`).
#' @param max_size_gb      Maximum size (GB) of the smoothed expression matrix
#'   that will be extracted to R (default `2.0`). Larger datasets print a
#'   warning and the matrix is left in Python memory only.
#' @param conda_env        Conda/virtualenv environment name.
#'
#' @return The Seurat object, with:
#'   - `Graphs[["esfs.knn"]]`: KNN distance graph (always).
#'   - Layer `"esfs.smoothed"` in the active assay: smoothed expression
#'     (only when the matrix is within `max_size_gb`).
#' @export
#' @seealso [GetESFSCellUMAPs()]
KNNSmoothESFS <- function(
  object,
  use_genes,
  assay            = Seurat::DefaultAssay(object),
  knn              = 30L,
  log_scale        = FALSE,
  force_recalculate = FALSE,
  max_size_gb      = 2.0,
  conda_env        = NULL
) {
  esfs  <- .get_esfs(conda_env)
  adata <- .build_and_restore_adata(object, assay = assay)

  adata <- esfs$knn_smooth_gene_expression(
    adata,
    use_genes         = reticulate::r_to_py(as.list(use_genes)),
    knn               = as.integer(knn),
    log_scale         = log_scale,
    force_recalculate = force_recalculate
  )

  # Always extract and store the KNN sparse graph (compact: n_cells × knn × 8 B)
  obsp_keys <- reticulate::py_to_r(reticulate::py_eval(
    "list(adata.obsp.keys())", convert = TRUE,
    local = list(adata = adata)
  ))
  if ("correlation_distance_kNN" %in% obsp_keys) {
    scipy   <- reticulate::import("scipy", convert = FALSE)
    knn_py  <- adata$obsp["correlation_distance_kNN"]
    knn_r   <- reticulate::py_to_r(knn_py)
    knn_graph <- as(knn_r, "Graph")
    rownames(knn_graph) <- colnames(object)
    colnames(knn_graph) <- colnames(object)
    object[["esfs.knn"]] <- knn_graph
  }

  # Extract smoothed expression only if within memory threshold
  n_cells  <- ncol(object)
  n_genes  <- nrow(object)
  size_gb  <- (n_cells * n_genes * 4) / 1e9

  if (size_gb <= max_size_gb) {
    smoothed <- reticulate::py_to_r(adata$layers["Smoothed_Expression"])
    # Python gives cells × genes; Seurat stores genes × cells
    smoothed_mat <- Matrix::Matrix(t(smoothed), sparse = FALSE)
    rownames(smoothed_mat) <- rownames(object)
    colnames(smoothed_mat) <- colnames(object)
    Seurat::LayerData(object, assay = assay, layer = "esfs.smoothed") <-
      smoothed_mat
    message(sprintf(
      "Smoothed expression (%.0f cells x %.0f genes, ~%.1f GB) stored in layer 'esfs.smoothed'.",
      n_cells, n_genes, size_gb
    ))
    message("Use FeaturePlot(object, features = \"GENE\", layer = \"esfs.smoothed\") to visualise.")
  } else {
    warning(sprintf(
      "Smoothed expression (~%.1f GB) exceeds max_size_gb=%.1f and was NOT transferred to R.\n",
      size_gb, max_size_gb
    ), "The matrix remains in Python memory (adata.layers['Smoothed_Expression']).\n",
    "Increase max_size_gb or access the result in Python directly.",
    call. = FALSE)
  }

  object <- Seurat::LogSeuratCommand(object)
  object
}


# ===========================================================================
# Visualisation
# ===========================================================================

#' Plot a UMAP of the top-ranked ESFS genes
#'
#' Wraps `esfs.plot_top_ranked_genes_UMAP()`. Embeds the top-ranked genes in
#' a 2-D UMAP based on their pairwise ESS correlation matrix and displays the
#' plot via matplotlib (rendered to the R graphics device as a PNG).
#'
#' @param object           A Seurat object with ES matrices and gene ranks.
#' @param top_ranked_genes Integer — how many top-ranked genes to embed.
#' @param assay            Assay to use.
#' @param label            Secondary features label (default `"Self"`).
#' @param clustering       Clustering method: `NULL` (no clustering),
#'   an integer (KMeans with that many clusters), or `"hdbscan"`.
#' @param known_genes      Character vector of gene names to highlight in the
#'   UMAP (plotted as black crosses).
#' @param umap_min_dist    UMAP `min_dist` parameter (default `0.1`).
#' @param umap_neighbours  UMAP `n_neighbors` parameter (default `20`).
#' @param hdbscan_min_size Minimum cluster size for HDBSCAN (default `50`).
#' @param random_state     Random seed for reproducibility (default `42L`).
#' @param dpi              Figure resolution in dots per inch (default `150L`).
#' @param conda_env        Conda/virtualenv environment name.
#'
#' @return Invisibly, a named list with:
#'   - `genes`: character vector of the top-ranked genes embedded
#'   - `labels`: integer cluster label per gene
#'   - `embedding`: numeric matrix (n_genes × 2) of UMAP coordinates
#' @export
#' @seealso [ESRankGenes()], [GetESFSCellUMAPs()]
PlotESFSGeneUMAP <- function(
  object,
  top_ranked_genes,
  assay            = Seurat::DefaultAssay(object),
  label            = "Self",
  clustering       = NULL,
  known_genes      = NULL,
  umap_min_dist    = 0.1,
  umap_neighbours  = 20L,
  hdbscan_min_size = 50L,
  random_state     = 42L,
  dpi              = 150L,
  conda_env        = NULL
) {
  esfs  <- .get_esfs(conda_env)
  np    <- reticulate::import("numpy", convert = FALSE)
  adata <- .build_and_restore_adata(object, assay = assay, label = label)

  py_clustering <- if (is.null(clustering) || identical(clustering, "None")) {
    reticulate::py_none()
  } else {
    clustering   # integer or "hdbscan" string
  }
  py_known <- if (is.null(known_genes)) {
    np$array(character(0))
  } else {
    np$array(known_genes)
  }

  result <- esfs$plot_top_ranked_genes_UMAP(
    adata,
    top_ranked_genes      = as.integer(top_ranked_genes),
    clustering            = py_clustering,
    known_important_genes = py_known,
    UMAP_min_dist         = umap_min_dist,
    UMAP_neighbours       = as.integer(umap_neighbours),
    hdbscan_min_cluster_size = as.integer(hdbscan_min_size),
    secondary_features_label = label,
    random_state          = as.integer(random_state)
  )

  .display_matplotlib_figure(dpi = dpi)

  out <- list(
    genes     = reticulate::py_to_r(result[[1]]),
    labels    = reticulate::py_to_r(result[[2]]),
    embedding = reticulate::py_to_r(result[[3]])
  )
  invisible(out)
}


#' Compute cell UMAP embeddings for each gene cluster
#'
#' Wraps `esfs.get_gene_cluster_cell_UMAPs()`. For each cluster of genes
#' identified by [PlotESFSGeneUMAP()], projects cells into a 2-D UMAP using
#' that cluster's gene set as features.
#'
#' **Backend notes (UMAP):**
#' - NVIDIA CUDA: KNN via cuML (O(N log N)) + cuML UMAP layout
#' - Apple Silicon / CPU: KNN via chunked BLAS matmul (O(N²)) + umap-learn
#'   (MLX does **not** accelerate UMAP)
#'
#' @param object              A Seurat object with gene rank and ES matrices.
#' @param gene_clust_labels   Integer vector of cluster labels, one per gene
#'   (as returned by [PlotESFSGeneUMAP()]).
#' @param top_ess_genes       Character vector of gene names corresponding to
#'   `gene_clust_labels` (as returned by [PlotESFSGeneUMAP()]).
#' @param assay               Assay to use.
#' @param n_neighbors         UMAP number of neighbours (default `30L`).
#' @param min_dist            UMAP `min_dist` (default `0.1`).
#' @param log_transformed     Log2-transform expression before UMAP
#'   (default `TRUE`).
#' @param specific_cluster    Integer or vector of integers — compute UMAP
#'   only for the specified cluster(s). `NULL` processes all clusters.
#' @param specific_genes      Character vector — bypass cluster-based
#'   selection and use these genes directly.
#' @param return_model        If `TRUE`, return fitted umap-learn models
#'   (required for [SaveESFSUMAPModel()]).
#' @param random_state        Random seed (default `NULL`).
#' @param conda_env           Conda/virtualenv environment name.
#'
#' @return A Python list `[embeddings, selected_genes]`.  Pass the result
#'   directly to [PlotESFSCellUMAPs()] or [SaveESFSUMAPModel()].
#' @export
#' @seealso [PlotESFSGeneUMAP()], [PlotESFSCellUMAPs()], [SaveESFSUMAPModel()]
GetESFSCellUMAPs <- function(
  object,
  gene_clust_labels,
  top_ess_genes,
  assay           = Seurat::DefaultAssay(object),
  n_neighbors     = 30L,
  min_dist        = 0.1,
  log_transformed = TRUE,
  specific_cluster = NULL,
  specific_genes  = NULL,
  return_model    = FALSE,
  random_state    = NULL,
  conda_env       = NULL
) {
  esfs  <- .get_esfs(conda_env)
  np    <- reticulate::import("numpy", convert = FALSE)
  adata <- .build_and_restore_adata(object, assay = assay)

  py_spec_cluster <- if (is.null(specific_cluster)) {
    reticulate::py_none()
  } else {
    reticulate::r_to_py(as.integer(specific_cluster))
  }
  py_spec_genes <- if (is.null(specific_genes)) {
    reticulate::py_none()
  } else {
    reticulate::r_to_py(as.list(specific_genes))
  }
  py_seed <- if (is.null(random_state)) {
    reticulate::py_none()
  } else {
    as.integer(random_state)
  }

  result <- esfs$get_gene_cluster_cell_UMAPs(
    adata,
    gene_clust_labels = np$array(as.integer(gene_clust_labels)),
    top_ESS_genes     = np$array(as.character(top_ess_genes)),
    n_neighbors       = as.integer(n_neighbors),
    min_dist          = min_dist,
    log_transformed   = log_transformed,
    specific_cluster  = py_spec_cluster,
    specific_genes    = py_spec_genes,
    return_model      = return_model,
    random_state      = py_seed
  )

  # Return Python objects directly — user passes them to PlotESFSCellUMAPs()
  # or SaveESFSUMAPModel() without R conversion.
  result
}


#' Plot cell UMAPs coloured by label or gene expression
#'
#' Wraps `esfs.plot_gene_cluster_cell_UMAPs()`. Takes the embeddings computed
#' by [GetESFSCellUMAPs()] and colours cells by a metadata label or gene
#' expression. All plotting is done in Python via matplotlib; the figure is
#' rendered to the R graphics device as a PNG.
#'
#' @param object                    A Seurat object.
#' @param gene_cluster_embeddings   Python list of embeddings (from
#'   [GetESFSCellUMAPs()]).
#' @param gene_cluster_selected_genes Python list of gene sets (from
#'   [GetESFSCellUMAPs()]).
#' @param assay                     Assay to use.
#' @param cell_label                Column in `object[[]]` (metadata) or gene
#'   name in `rownames(object)` used to colour cells. Use `"None"` for no
#'   colouring.
#' @param ncol                      Number of plot columns (default `1`).
#' @param log2_gene_expression      Log2-transform expression when colouring
#'   by a gene (default `TRUE`).
#' @param figsize                   Tuple `c(width, height)` in inches
#'   (default `c(18, 10)`).
#' @param marker_size               Point size (default `3L`).
#' @param sort_by_value             Sort cells by expression before plotting
#'   so high-expressing cells appear on top (default `TRUE`).
#' @param dpi                       Figure resolution (default `150L`).
#' @param conda_env                 Conda/virtualenv environment name.
#'
#' @return Invisibly `NULL`. The figure is rendered to the active R graphics device.
#' @export
#' @seealso [GetESFSCellUMAPs()]
PlotESFSCellUMAPs <- function(
  object,
  gene_cluster_embeddings,
  gene_cluster_selected_genes,
  assay               = Seurat::DefaultAssay(object),
  cell_label          = "None",
  ncol                = 1L,
  log2_gene_expression = TRUE,
  figsize             = c(18, 10),
  marker_size         = 3L,
  sort_by_value       = TRUE,
  dpi                 = 150L,
  conda_env           = NULL
) {
  esfs  <- .get_esfs(conda_env)
  adata <- .build_and_restore_adata(object, assay = assay)

  esfs$plot_gene_cluster_cell_UMAPs(
    adata,
    gene_cluster_embeddings    = gene_cluster_embeddings,
    gene_cluster_selected_genes = gene_cluster_selected_genes,
    cell_label           = cell_label,
    ncol                 = as.integer(ncol),
    log2_gene_expression = log2_gene_expression,
    figsize              = reticulate::tuple(as.numeric(figsize)),
    marker_size          = as.integer(marker_size),
    sort_by_value        = sort_by_value
  )

  .display_matplotlib_figure(dpi = dpi)
  invisible(NULL)
}


# ===========================================================================
# UMAP model persistence
# ===========================================================================

#' Save a fitted UMAP model to disk
#'
#' Wraps `esfs.save_umap_model()`. Saves the model, gene list, and metadata
#' to a single `.joblib` file. The saved model can be used to project new
#' cells into the same embedding with [LoadESFSUMAPModel()].
#'
#' @param model           A fitted UMAP model (Python object from
#'   [GetESFSCellUMAPs()] with `return_model = TRUE`).
#' @param gene_list       Character vector of genes used to fit the model.
#' @param filepath        Output file path (`.joblib` extension recommended).
#' @param log_transformed Whether the training data was log2(X+1) transformed
#'   (default `FALSE`).
#' @param conda_env       Conda/virtualenv environment name.
#'
#' @return Invisibly `NULL`.
#' @export
#' @seealso [LoadESFSUMAPModel()], [GetESFSCellUMAPs()]
SaveESFSUMAPModel <- function(
  model,
  gene_list,
  filepath,
  log_transformed = FALSE,
  conda_env       = NULL
) {
  esfs <- .get_esfs(conda_env)
  np   <- reticulate::import("numpy", convert = FALSE)

  esfs$save_umap_model(
    model,
    np$array(as.character(gene_list)),
    filepath,
    log_transformed
  )
  invisible(NULL)
}


#' Load a saved UMAP model
#'
#' Wraps `esfs.load_umap_model()`. Loads a model bundle previously saved with
#' [SaveESFSUMAPModel()].
#'
#' @param filepath  Path to the `.joblib` file.
#' @param conda_env Conda/virtualenv environment name.
#'
#' @return A Python dict with keys:
#'   - `model`: the fitted umap-learn UMAP model (supports `.transform()`)
#'   - `gene_list`: genes used for the embedding
#'   - `log_transformed`: whether log2(X+1) was applied before fitting
#' @export
#' @seealso [SaveESFSUMAPModel()]
LoadESFSUMAPModel <- function(filepath, conda_env = NULL) {
  esfs <- .get_esfs(conda_env)
  esfs$load_umap_model(filepath)
}


# ===========================================================================
# Backend configuration
# ===========================================================================

#' Configure the ESFS computation backend
#'
#' Wraps `esfs.configure()`. Sets the hardware backend (CPU, CUDA GPU, or
#' Apple Silicon MLX) and floating-point precision for the ES matrix
#' calculations.
#'
#' @param gpu     If `TRUE`, attempt to use GPU (CUDA first, then MLX).
#' @param upcast  If `TRUE`, use float64 precision (uses more memory).
#'   Default `FALSE` (float32).
#' @param verbose If `TRUE`, print the backend status banner to the console.
#' @param conda_env Conda/virtualenv environment name.
#'
#' @return Invisibly `NULL`.
#' @export
#' @seealso [UseESFSGPU()], [UseESFSMLX()], [UseESFSCPU()], [GetESFSBackendInfo()]
ConfigureESFSBackend <- function(
  gpu       = TRUE,
  upcast    = FALSE,
  verbose   = TRUE,
  conda_env = NULL
) {
  esfs <- .get_esfs(conda_env)
  esfs$configure(gpu = gpu, upcast = upcast, verbose = verbose)
  invisible(NULL)
}


#' Force ESFS to use the CPU backend
#'
#' @param conda_env Conda/virtualenv environment name.
#' @return Invisibly `NULL`.
#' @export
#' @seealso [ConfigureESFSBackend()], [UseESFSGPU()], [UseESFSMLX()]
UseESFSCPU <- function(conda_env = NULL) {
  esfs <- .get_esfs(conda_env)
  esfs$use_cpu()
  invisible(NULL)
}


#' Force ESFS to use the GPU backend (CUDA or MLX)
#'
#' Tries CUDA (CuPy) first, then MLX (Apple Silicon), then falls back to CPU
#' with a warning if neither is available.
#'
#' @param conda_env Conda/virtualenv environment name.
#' @return Invisibly `NULL`.
#' @export
#' @seealso [ConfigureESFSBackend()], [UseESFSCPU()], [UseESFSMLX()]
UseESFSGPU <- function(conda_env = NULL) {
  esfs <- .get_esfs(conda_env)
  esfs$use_gpu()
  invisible(NULL)
}


#' Force ESFS to use the Apple Silicon MLX/Metal backend
#'
#' Only effective on Apple Silicon Macs (M1 and later) with MLX installed.
#' Falls back to CPU with a warning on other hardware.
#'
#' @param conda_env Conda/virtualenv environment name.
#' @return Invisibly `NULL`.
#' @export
#' @seealso [ConfigureESFSBackend()], [UseESFSCPU()], [UseESFSGPU()]
UseESFSMLX <- function(conda_env = NULL) {
  esfs <- .get_esfs(conda_env)
  esfs$use_mlx()
  invisible(NULL)
}


#' Report the active ESFS computation backend
#'
#' @param conda_env Conda/virtualenv environment name.
#' @return A character string describing the active backend and precision.
#' @export
#' @seealso [ConfigureESFSBackend()]
GetESFSBackendInfo <- function(conda_env = NULL) {
  esfs <- .get_esfs(conda_env)
  reticulate::py_to_r(esfs$get_backend_info())
}
