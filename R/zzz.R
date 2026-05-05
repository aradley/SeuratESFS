# Package-level environment for caching Python module references
.esfs_pkg_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Resolve the configured environment name
  env <- getOption("esfs.conda_env",
                   Sys.getenv("ESFS_CONDA_ENV", unset = NA_character_))
  if (is.na(env)) env <- "r-esfs"

  # Only configure reticulate if Python has not yet been initialised.
  # py_available(initialize = FALSE) peeks without starting the interpreter.
  if (!reticulate::py_available(initialize = FALSE)) {
    if (reticulate::condaenv_exists(env)) {
      reticulate::use_condaenv(env, required = FALSE)
    } else if (reticulate::virtualenv_exists(env)) {
      reticulate::use_virtualenv(env, required = FALSE)
    }
    # If neither exists, stay silent here — .get_esfs() will surface a clear
    # error with fix instructions the first time a wrapper function is called.
  }
}

.onAttach <- function(libname, pkgname) {
  env <- getOption("esfs.conda_env",
                   Sys.getenv("ESFS_CONDA_ENV", unset = "r-esfs"))

  env_found <- reticulate::condaenv_exists(env) ||
               reticulate::virtualenv_exists(env)

  if (!env_found) {
    packageStartupMessage(
      "SeuratESFS: Python environment '", env, "' not found.\n",
      "Run install_esfs() to create it automatically.\n",
      "If esfs is already installed elsewhere, set before loading this package:\n",
      "  options(esfs.conda_env = \"your-env-name\")"
    )
  }
}
