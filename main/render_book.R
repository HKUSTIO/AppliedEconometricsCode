bookdown::render_book(
  input = "bookdown",
  output_format = "bookdown::gitbook",
  preview = TRUE
)

bookdown::serve_book(
  dir = "bookdown",
  preview = TRUE
)

bookdown::render_book(
  input = "bookdown",
  output_format = "bookdown::pdf_book",
  preview = TRUE
)
