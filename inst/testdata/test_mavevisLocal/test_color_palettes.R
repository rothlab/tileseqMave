#!/usr/bin/env Rscript

# Smoke tests for genophenogram color palette handling.
# This script shells out to the actual CLI wrapper so it exercises the same
# argument parsing and plotting path as a real user run.

options(stringsAsFactors = FALSE)

resolve_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(dirname(sys.frames()[[1]]$ofile), winslash = "/", mustWork = TRUE))
  }

  stop("Could not determine the script directory. Run this script directly with Rscript.")
}

script_dir <- resolve_script_dir()
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = TRUE)
workspace <- script_dir
input_dir <- file.path(workspace, "demo_scores")
param_file <- file.path(workspace, "parameters.json")
out_dir <- file.path(workspace, "output")
wrapper <- file.path(repo_root, "inst/scripts/mavevisLocal.R")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_case <- function(label, palette = NULL, custom_colors = NULL) {
  out_file <- file.path(out_dir, paste0("test_", label, ".pdf"))

  cmd <- paste(
    "Rscript",
    shQuote(wrapper),
    "--workspace", shQuote(workspace),
    "--input", shQuote(input_dir),
    "--parameters", shQuote(param_file),
    "--output", shQuote(out_file)
  )

  if (!is.null(palette)) {
    cmd <- paste(cmd, "--colorPalette", shQuote(palette))
  }
  if (!is.null(custom_colors)) {
    cmd <- paste(cmd, "--customColors", shQuote(custom_colors))
  }

  cat("\n=== Running:", label, "===\n")
  cat(cmd, "\n")

  status <- system(cmd)
  if (!file.exists(out_file)) {
    stop("Output PDF was not created for: ", label)
  }
  if (status != 0) {
    stop("Palette smoke test failed for: ", label)
  }

  cat("OK:", label, "->", out_file, "\n")
}

run_case("default", palette = "default")
run_case("viridis", palette = "viridis")
run_case("custom", palette = "custom", custom_colors = "#1B2A41,#D9D9D9,#D1495B")

# Invalid palette should fail
cat("\n=== Running: invalid_palette_should_fail ===\n")
invalid_file <- file.path(out_dir, "test_invalid_palette.pdf")
cmd <- paste(
  "Rscript",
  shQuote(wrapper),
  "--workspace", shQuote(workspace),
  "--input", shQuote(input_dir),
  "--parameters", shQuote(param_file),
  "--output", shQuote(invalid_file),
  "--colorPalette", shQuote("notARealPalette")
)

status <- system(cmd)
if (status == 0) {
  stop("Invalid palette test did not fail as expected.")
}

cat("\nAll palette smoke tests passed.\n")
