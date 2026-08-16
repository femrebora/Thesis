# .Rprofile — Thesis reproducible analysis repository
# Loads renv only when a lockfile already exists.

if (file.exists("renv.lock")) {
  source("renv/activate.R")
} else {
  message(
    "No renv.lock present. ENVIRONMENT.txt is the historical environment ",
    "manifest for the thesis figure run. Do not run renv::init() here to ",
    "\"capture\" the environment — that would record today's packages, not ",
    "the thesis versions. Verify with: Rscript tools/check_environment.R"
  )
}
