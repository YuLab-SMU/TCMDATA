#!/usr/bin/env Rscript
# Build script for bookdown - used by GitHub Actions CI/CD

# Set locale for CJK character support
tryCatch(
  Sys.setlocale("LC_ALL", "en_US.UTF-8"),
  warning = function(w) message("Locale warning: ", w$message),
  error   = function(e) message("Locale error: ", e$message)
)

message("=== Build Environment ===")
message("R version: ", R.version.string)
message("Working dir: ", getwd())
message("Locale: ", Sys.getlocale())
message("Files in docs/: ", paste(list.files(), collapse = ", "))

# Check bookdown is available
if (!requireNamespace("bookdown", quietly = TRUE)) {
  stop("bookdown package not found!")
}
message("bookdown version: ", packageVersion("bookdown"))

# Set output directory
output_dir <- file.path(normalizePath(".."), "site")
message("Output dir: ", output_dir)

# Render the book
message("=== Starting bookdown render ===")
bookdown::render_book(
  input = "index.Rmd",
  output_format = "bookdown::gitbook",
  output_dir = output_dir
)
message("=== Build complete ===")
