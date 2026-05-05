#' Install the ESFS Python environment
#'
#' Creates a dedicated Python environment and installs the ESFS package along
#' with all required dependencies. Safe to call repeatedly:
#'
#' - **First run**: creates the `r-esfs` conda environment and installs all
#'   dependencies (~5 minutes).
#' - **Subsequent runs** (e.g. after updating SeuratESFS via
#'   `remotes::install_github("aradley/SeuratESFS")`): detects the existing
#'   environment and only force-reinstalls the ESFS Python package from GitHub
#'   (~30 seconds), leaving all other dependencies untouched.
#'
#' The `~/.Rprofile` option `esfs.conda_env` is written automatically so the
#' correct Python environment is found in every future R session.
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

  extras   <- switch(backend, cpu = "", mlx = "[mlx]", gpu = "[gpu]")
  esfs_pkg <- paste0(
    "esfs", extras,
    " @ git+https://github.com/aradley/ESFS.git@memory_optimised"
  )

  # ---- Fast path: environment already exists --------------------------------
  # Just force-reinstall ESFS from GitHub; skip recreating the env and
  # reinstalling all core deps (saves ~5 minutes on subsequent calls).
  conda_exists <- reticulate::condaenv_exists(envname)
  venv_exists  <- reticulate::virtualenv_exists(envname)

  if (conda_exists || venv_exists) {
    existing_method <- if (conda_exists) "conda" else "virtualenv"
    message("Environment '", envname, "' already exists.")
    message("Reinstalling ESFS from GitHub (memory_optimised branch)...")
    message("  Backend : ", backend)
    reticulate::py_install(
      esfs_pkg,
      envname = envname,
      method  = existing_method,
      pip     = TRUE
    )
    .verify_esfs_installation(envname, existing_method)
    message("\n=== Reinstall complete ===")
    .add_rprofile_option(envname)
    message("Restart R to pick up the updated ESFS package.")
    if (restart_session &&
          requireNamespace("rstudioapi", quietly = TRUE) &&
          rstudioapi::isAvailable()) {
      rstudioapi::restartSession()
    }
    return(invisible(envname))
  }

  # ---- Full path: create environment from scratch --------------------------
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

  if (method == "conda") {
    reticulate::conda_create(envname        = envname,
                             python_version = python_version)
  } else {
    reticulate::virtualenv_create(envname = envname,
                                  python  = python_version)
  }

  # numba/llvmlite must be installed via conda, not pip.
  # pip-installed llvmlite links against libc++.1.dylib via @rpath, which
  # macOS cannot resolve when Python is loaded through R/reticulate because
  # R.framework ends up in the rpath and libc++.1.dylib is not there.
  # conda install pulls libcxx as a proper conda dep, placing
  # libc++.1.dylib in the env's lib/ where the linker finds it.
  if (method == "conda") {
    message(
      "\nInstalling numba via conda (required for macOS + R compatibility)..."
    )
    reticulate::conda_install(
      envname  = envname,
      packages = "numba>=0.60,<0.62",
      channel  = "conda-forge",
      pip      = FALSE
    )
  }

  # Everything else can be pip-installed (numba excluded for conda path)
  core_pkgs <- c(
    "numpy", "scipy<2.0", "pandas", "matplotlib", "plotly",
    "anndata>=0.10,<0.12", "scanpy[leiden]>=1.10,<2.0",
    "scikit-learn>=1.4,<2.0", "umap-learn>=0.5,<0.6",
    "multiprocess>=0.70,<0.80",
    "p-tqdm>=1.4,<2.0", "tqdm", "zarr>=2.16,<3.0"
  )
  if (method != "conda") {
    core_pkgs <- c(core_pkgs, "numba>=0.60,<0.62")
  }

  message("\nInstalling core dependencies...")
  reticulate::py_install(core_pkgs, envname = envname,
                         method = method, pip = TRUE)

  message(
    "\nInstalling ESFS", extras,
    " from GitHub (memory_optimised branch)..."
  )
  reticulate::py_install(esfs_pkg, envname = envname,
                         method = method, pip = TRUE)

  .verify_esfs_installation(envname, method)

  message("\n=== Installation complete ===")
  .add_rprofile_option(envname)
  message("Restart R.")

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
  seurat_esfs_ver <- tryCatch(
    as.character(utils::packageVersion("SeuratESFS")),
    error = function(e) "?"
  )
  cat("[R] SeuratESFS :", seurat_esfs_ver, "\n")
  cat("[R] reticulate :",
      as.character(utils::packageVersion("reticulate")), "\n\n")

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
    cfg     <- reticulate::py_config()
    # cfg$version can be a list in newer reticulate — flatten safely
    ver_str <- tryCatch(
      paste(unlist(cfg[["version"]]), collapse = "."),
      error = function(e) "unknown"
    )
    cat("[OK] Python ", ver_str, "\n", sep = "")
    cat("     Path   : ", cfg$python, "\n", sep = "")

    # Warn if the active Python is not inside the expected environment
    expected_env_path <- tryCatch({
      if (reticulate::condaenv_exists(envname))
        reticulate::conda_python(envname = envname)
      else
        file.path(reticulate::virtualenv_root(), envname, "bin", "python")
    }, error = function(e) NULL)

    active   <- normalizePath(cfg$python, mustWork = FALSE)
    expected <- normalizePath(expected_env_path %||% "", mustWork = FALSE)
    if (!is.null(expected_env_path) &&
          !identical(active, expected)) {
      cat("[!!] WARNING: Active Python is NOT from '", envname,
          "'!\n", sep = "")
      cat("     Expected: ", expected_env_path, "\n", sep = "")
      cat("     Fix: add  options(esfs.conda_env = \"", envname,
          "\")  to ~/.Rprofile\n", sep = "")
      cat("          then restart R and load SeuratESFS BEFORE any other\n")
      cat("          package that uses Python.\n")
    }
    cat("\n")
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
.add_rprofile_option <- function(envname) {
  rprofile    <- path.expand("~/.Rprofile")
  option_line <- paste0('options(esfs.conda_env = "', envname, '")')

  existing <- if (file.exists(rprofile)) {
    readLines(rprofile, warn = FALSE)
  } else {
    character(0)
  }

  if (any(grepl("esfs.conda_env", existing, fixed = TRUE))) {
    message(
      "  ~/.Rprofile already contains esfs.conda_env — no change needed."
    )
    return(invisible(NULL))
  }

  # Append with a preceding blank line for readability
  write(c("", option_line), file = rprofile, append = TRUE)
  message("  Added to ~/.Rprofile: ", option_line)
  message(
    "  This ensures the correct Python environment is used in every session."
  )
}

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

  # If Python is already running we cannot switch environments.
  if (reticulate::py_available(initialize = FALSE)) {
    message(
      "Note: Python is already initialised in this session; skipping\n",
      "  environment switch. Run check_esfs() in a fresh R session to confirm."
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
  }, error = function(e) {
    msg <- conditionMessage(e)
    message("[!!] esfs import failed: ", msg)
    if (grepl("shared object|dylib|OSError|LoadLibrary", msg,
              ignore.case = TRUE)) {
      message(
        "     Likely cause: the GPU/Metal backend library failed to load.\n",
        "     Fix options:\n",
        "       1. CPU-only (recommended): delete the env and reinstall:\n",
        "            conda env remove -n ", envname,
        "  # run in terminal\n",
        "          then: install_esfs(backend = \"cpu\")\n",
        "       2. Diagnose: run  reticulate::py_last_error()  for details."
      )
    } else {
      message(
        "     Run check_esfs() or reticulate::py_last_error() for details."
      )
    }
  })
}
