#' Install the ESFS Python environment
#'
#' Creates a dedicated Python environment and installs the ESFS package along
#' with all required dependencies. This function needs to be called only once
#' per machine. After installation, add the following line to your
#' `~/.Rprofile` so the environment is found automatically in every session:
#'
#' ```r
#' options(esfs.conda_env = "r-esfs")
#' ```
#'
#' @param envname Name of the conda or virtualenv environment to create.
#'   Defaults to `"r-esfs"`.
#' @param backend Hardware acceleration backend. One of `"auto"` (detect from
#'   hardware), `"cpu"` (NumPy/Numba only), `"mlx"` (Apple Silicon Metal),
#'   or `"gpu"` (NVIDIA CUDA via CuPy).
#' @param python_version Python version to install. Defaults to `"3.11"`, which
#'   is compatible with all optional backends including cuML.
#' @param method Environment manager. One of `"conda"` (preferred — more
#'   reliable for numba) or `"virtualenv"`. Falls back to virtualenv if conda
#'   is not available.
#' @param restart_session If `TRUE` and running inside RStudio, restart the
#'   R session after installation so the new environment is picked up cleanly.
#'
#' @return Invisibly returns the environment name.
#' @export
#'
#' @examples
#' \dontrun{
#' # First-time setup — auto-detects MLX on Apple Silicon, CUDA on NVIDIA
#' install_esfs()
#'
#' # Force CPU-only installation
#' install_esfs(backend = "cpu")
#'
#' # Custom environment name
#' install_esfs(envname = "my-esfs-env")
#' }
install_esfs <- function(
  envname         = "r-esfs",
  backend         = c("auto", "cpu", "mlx", "gpu"),
  python_version  = "3.11",
  method          = c("conda", "virtualenv"),
  restart_session = TRUE
) {
  backend <- match.arg(backend)
  method  <- match.arg(method)

  if (backend == "auto") backend <- .detect_esfs_backend()

  message("Installing ESFS Python environment")
  message("  Environment : ", envname)
  message("  Backend     : ", backend)
  message("  Python      : ", python_version)
  message("  Method      : ", method)

  # Fall back to virtualenv if conda is unavailable
  if (method == "conda") {
    conda_bin <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
    if (is.null(conda_bin) || !nzchar(conda_bin)) {
      message("Conda not found — falling back to virtualenv.")
      method <- "virtualenv"
    }
  }

  # Create the isolated environment
  if (method == "conda") {
    reticulate::conda_create(envname        = envname,
                             python_version = python_version)
  } else {
    reticulate::virtualenv_create(envname = envname,
                                  python  = python_version)
  }

  # Core Python packages (always required)
  core_pkgs <- c(
    "numpy", "scipy<2.0", "pandas", "matplotlib", "plotly",
    "anndata>=0.10,<0.12", "scanpy[leiden]>=1.10,<2.0",
    "scikit-learn>=1.4,<2.0", "umap-learn>=0.5,<0.6",
    "numba>=0.60,<0.62", "multiprocess>=0.70,<0.80",
    "p-tqdm>=1.4,<2.0", "tqdm", "zarr>=2.16,<3.0"
  )

  message("\nInstalling core dependencies...")
  reticulate::py_install(core_pkgs, envname = envname,
                         method = method, pip = TRUE)

  # ESFS package from GitHub with optional extras
  extras   <- switch(backend, cpu = "", mlx = "[mlx]", gpu = "[gpu]")
  esfs_pkg <- paste0(
    "esfs", extras,
    " @ git+https://github.com/aradley/ESFS.git@memory_optimised"
  )

  message("\nInstalling ESFS", extras, " from GitHub...")
  reticulate::py_install(esfs_pkg, envname = envname,
                         method = method, pip = TRUE)

  # Verify everything works
  .verify_esfs_installation(envname, method)

  message("\n=== Installation complete ===")
  message("Add the following line to your ~/.Rprofile for automatic setup:")
  message('  options(esfs.conda_env = "', envname, '")')
  message("Then restart R.")

  if (restart_session &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    rstudioapi::restartSession()
  }

  invisible(envname)
}


#' Check the SeuratESFS installation
#'
#' Runs a diagnostic report showing the status of the Python environment,
#' installed packages, and active ESFS backend. Run this first when
#' troubleshooting installation issues.
#'
#' @param envname Environment name to check. Defaults to the value of
#'   `getOption("esfs.conda_env")` or `"r-esfs"`.
#'
#' @return Invisibly returns `NULL`.
#' @export
#'
#' @examples
#' \dontrun{
#' check_esfs()
#' }
check_esfs <- function(envname = NULL) {
  envname <- envname %||%
    getOption("esfs.conda_env",
              Sys.getenv("ESFS_CONDA_ENV", unset = "r-esfs"))

  cat("=== SeuratESFS Diagnostic Report ===\n\n")

  # R package versions
  cat("[R] SeuratESFS :", tryCatch(
        as.character(utils::packageVersion("SeuratESFS")), error = function(e) "?"),
      "\n")
  cat("[R] reticulate :", as.character(utils::packageVersion("reticulate")), "\n\n")

  # Environment existence
  conda_ok <- reticulate::condaenv_exists(envname)
  venv_ok  <- reticulate::virtualenv_exists(envname)
  if (conda_ok) {
    cat("[OK] Conda env '", envname, "' found\n", sep = "")
  } else if (venv_ok) {
    cat("[OK] Virtualenv '", envname, "' found\n", sep = "")
  } else {
    cat("[!!] Environment '", envname, "' NOT found\n", sep = "")
    cat("     Run install_esfs() to create it.\n\n")
    return(invisible(NULL))
  }

  # Python interpreter
  if (reticulate::py_available()) {
    cfg <- reticulate::py_config()
    cat("[OK] Python ", cfg$version, "\n", sep = "")
    cat("     Path: ", cfg$python, "\n\n", sep = "")
  } else {
    cat("[!!] Python not initialised\n\n")
    return(invisible(NULL))
  }

  # Core Python packages
  cat("[Core packages]\n")
  core_mods <- c("esfs", "anndata", "numpy", "scipy", "numba",
                 "umap", "sklearn", "matplotlib", "scanpy")
  for (pkg in core_mods) {
    ok <- reticulate::py_module_available(pkg)
    cat(if (ok) "[OK] " else "[!!] ", pkg, "\n", sep = "")
  }

  # Optional backends
  cat("\n[Optional backends]\n")
  for (pkg in c("cupy", "mlx")) {
    ok <- reticulate::py_module_available(pkg)
    cat(if (ok) "[OK] " else "[ ] ", pkg,
        if (!ok) "  (optional — not installed)" else "",
        "\n", sep = "")
  }

  # Active ESFS backend
  if (reticulate::py_module_available("esfs")) {
    esfs <- reticulate::import("esfs", convert = FALSE)
    esfs$configure(gpu = TRUE, upcast = FALSE, verbose = FALSE)
    info <- reticulate::py_to_r(esfs$get_backend_info())
    cat("\n[ESFS] ", info, "\n", sep = "")
  }

  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
.detect_esfs_backend <- function() {
  sysinfo <- Sys.info()

  # Apple Silicon
  if (identical(sysinfo[["sysname"]], "Darwin")) {
    machine <- tryCatch(trimws(system("uname -m", intern = TRUE)),
                        error = function(e) "")
    if (identical(machine, "arm64")) {
      message("  Detected: Apple Silicon — MLX/Metal backend")
      return("mlx")
    }
  }

  # NVIDIA GPU
  if (nzchar(Sys.which("nvidia-smi"))) {
    message("  Detected: NVIDIA GPU — CUDA backend")
    return("gpu")
  }

  message("  Detected: CPU only")
  "cpu"
}

#' @keywords internal
.verify_esfs_installation <- function(envname, method = "conda") {
  message("\nVerifying installation...")

  # If Python is already running we cannot switch environments — skip the env
  # switch and advise the user to verify in a fresh session instead.
  if (reticulate::py_available(initialize = FALSE)) {
    message(
      "Note: Python is already initialised in this session; skipping\n",
      "  environment switch. The new '", envname, "' environment will be\n",
      "  active after you restart R. Run check_esfs() to confirm."
    )
    return(invisible(NULL))
  }

  if (method == "conda") {
    reticulate::use_condaenv(envname, required = TRUE)
  } else {
    reticulate::use_virtualenv(envname, required = TRUE)
  }

  tryCatch({
    esfs <- reticulate::import("esfs")
    ver  <- reticulate::py_to_r(esfs$`__version__`)
    message("[OK] esfs ", ver, " imported successfully")
    esfs$configure(gpu = TRUE, upcast = FALSE, verbose = FALSE)
    message("[OK] Backend: ", reticulate::py_to_r(esfs$get_backend_info()))
  }, error = function(e) {
    stop(
      "Installation verification failed: ", conditionMessage(e), "\n",
      "Run check_esfs() for a full diagnostic report."
    )
  })
}
