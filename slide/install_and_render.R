user_lib <- Sys.getenv("R_LIBS_USER")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

pkgs <- c("rmarkdown", "knitr", "magrittr", "purrr", "foreach", "ggplot2", "modelsummary", "kableExtra", "haven", "rdrobust", "rdpower")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org", lib = user_lib)
  }
}

setwd(dirname(normalizePath("slide")))
files_rmd <- list.files(path = "slide", pattern = "\\.Rmd$", full.names = TRUE)
for (f in files_rmd) {
  message("Rendering ", f)
  result <- try(rmarkdown::render(f), silent = FALSE)
  if (inherits(result, "try-error")) {
    warning("Failed to render ", f, ": ", attr(result, "condition")$message)
  }
}
