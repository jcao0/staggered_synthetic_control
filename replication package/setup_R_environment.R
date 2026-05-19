## Optional R setup for exact replication of the archived GSC results.
##
## The GSC results in Figure 3 were generated with gsynth 1.2.1.
## Newer versions of gsynth changed implementation/default behavior and can
## produce different GSC estimates, especially for Figure 3 panels e and f.
##
## This script attempts to install the archived version used for replication.
## If installation is not possible on a user's machine, the replication scripts
## can still be run with a newer gsynth version, but the GSC results may differ.

required_gsynth_version <- "1.2.1"

options(repos = c(CRAN = "https://cran.r-project.org"))

install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

message("This setup script attempts to install gsynth ", required_gsynth_version,
        " for exact replication of the archived GSC results.")
message("Users may skip this script and use a newer gsynth version, but Figure 3 GSC results may differ.")

tryCatch(
  install_if_missing("remotes"),
  error = function(e) {
    warning(
      "Could not install remotes: ", conditionMessage(e),
      "\nThe script will continue without installing archived gsynth.",
      call. = FALSE
    )
  }
)

current_gsynth_version <- if (requireNamespace("gsynth", quietly = TRUE)) {
  as.character(utils::packageVersion("gsynth"))
} else {
  NA_character_
}

if (!identical(current_gsynth_version, required_gsynth_version) &&
    requireNamespace("remotes", quietly = TRUE)) {
  tryCatch(
    {
      remotes::install_version(
        "gsynth",
        version = required_gsynth_version,
        repos = "https://cran.r-project.org",
        upgrade = "never",
        dependencies = TRUE
      )
    },
    error = function(e) {
      warning(
        "Could not install gsynth ", required_gsynth_version, ": ",
        conditionMessage(e),
        "\nThe replication scripts can still run with the installed gsynth version, ",
        "but Figure 3 GSC results may differ from the archived outputs.",
        call. = FALSE
      )
    }
  )
} else if (!requireNamespace("remotes", quietly = TRUE)) {
  warning(
    "Package remotes is not available, so gsynth ",
    required_gsynth_version,
    " was not installed automatically. The replication scripts can still run ",
    "with the installed gsynth version, but Figure 3 GSC results may differ.",
    call. = FALSE
  )
}

for (package in c("dplyr", "panelView", "ggplot2", "cowplot")) {
  install_if_missing(package)
}

if (requireNamespace("gsynth", quietly = TRUE)) {
  message("Installed gsynth version: ", as.character(utils::packageVersion("gsynth")))
} else {
  warning("gsynth is not installed. Please install gsynth before running estimation_gsc.R.", call. = FALSE)
}
