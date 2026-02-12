
# initialize --------------------------------------------------------------

library(magrittr)
library(foreach)



# read files --------------------------------------------------------------

files_rmd <-
  list.files(
    path = "slide",
    pattern = "Rmd",
    full.names = TRUE
  )


# knit files --------------------------------------------------------------

# inside the project folder
files_rmd %>%
  purrr::map(
    rmarkdown::render
  )

# outside the project folder
dir.create("../AppliedEconometricsPublic", showWarnings = FALSE)
for (i in 1:length(files_rmd)) {
    rmarkdown::render(
      files_rmd[i],
      output_dir = "../AppliedEconometricsPublic"
    )
}
