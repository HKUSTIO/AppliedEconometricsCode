file_list <-
  list.files(
    path = "main",
    pattern = "Rmd",
    full.names = TRUE
  )

for (
  filename in file_list
  ) {
    savename <-
      gsub(
        "Rmd",
        "html",
        filename
      )
    rmarkdown::render(
      input = 
        filename %>% 
        here::here(),
      output_file = 
        savename %>% 
        here::here()
    )

  }
