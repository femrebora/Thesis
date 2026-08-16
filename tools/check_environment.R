#!/usr/bin/env Rscript
# =============================================================================
# check_environment.R — compare the installed R packages against ENVIRONMENT.txt
# =============================================================================
# The thesis figures were produced with the exact versions recorded in
# ENVIRONMENT.txt. Package updates — particularly ggplot2 and scales, which own
# axis-break selection and panel layout — can shift a figure without changing
# any underlying value. Run this before trusting a regenerated figure set.
#
#   Rscript tools/check_environment.R
#
# Exit status 0 = R version and every listed package match.
# Exit status 1 = R version differs, and/or at least one package differs or is
# missing.
# =============================================================================

manifest <- here::here("ENVIRONMENT.txt")
if (!file.exists(manifest)) {
  stop("ENVIRONMENT.txt not found at ", manifest, call. = FALSE)
}

lines <- readLines(manifest, warn = FALSE)

r_recorded <- sub("^R version: *", "", grep("^R version:", lines, value = TRUE)[1])
r_current  <- R.version.string
r_ok <- identical(r_recorded, r_current)
cat(sprintf("%-22s %s\n", "R", if (r_ok) "match"
            else paste0("DIFFERS\n    recorded: ", r_recorded,
                        "\n    current : ", r_current)))

# Package rows follow the PACKAGE/VERSION header line.
hdr <- grep("^PACKAGE", lines)
rows <- if (length(hdr)) lines[(hdr[1] + 1):length(lines)] else character()
rows <- rows[nzchar(trimws(rows))]

mismatches <- 0L
for (row in rows) {
  parts <- strsplit(trimws(row), " +")[[1]]
  if (length(parts) < 2) next
  pkg <- parts[1]
  recorded <- parts[2]
  current <- tryCatch(as.character(utils::packageVersion(pkg)),
                      error = function(e) NA_character_)

  if (is.na(current)) {
    cat(sprintf("%-22s NOT INSTALLED (recorded %s)\n", pkg, recorded))
    mismatches <- mismatches + 1L
  } else if (!identical(current, recorded)) {
    cat(sprintf("%-22s DIFFERS  recorded %s -> installed %s\n",
                pkg, recorded, current))
    mismatches <- mismatches + 1L
  } else {
    cat(sprintf("%-22s %s\n", pkg, current))
  }
}

failed <- (!r_ok) || (mismatches > 0L)

if (mismatches > 0L) {
  cat("\n", mismatches, " package(s) differ from the recorded environment.\n", sep = "")
  cat("Figures may still be scientifically correct, but can differ cosmetically\n")
  cat("(axis breaks, panel layout, label placement). Compare PNGs before use:\n")
  cat("  git status --short -- '*.png'\n")
}
if (!r_ok) {
  cat("\nR version differs from the recorded thesis environment.\n")
  cat("  recorded: ", r_recorded, "\n", sep = "")
  cat("  current : ", r_current, "\n", sep = "")
}
if (!failed) {
  cat("\nEnvironment matches the one that produced the thesis figures.\n")
}

quit(save = "no", status = if (failed) 1L else 0L)
