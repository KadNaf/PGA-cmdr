# ── 01_start.R ─────────────────────────────────────────────────────────────
# Run once to initialise the golem package.
# Do NOT source this file as part of the app.

# 1) Install / update golem if needed
# install.packages("golem")

# 2) Set active project to pgacmdr/
# usethis::proj_set(".")

# 3) Fill in DESCRIPTION fields (run interactively)
# golem::fill_desc(
#   pkg_name  = "pgacmdr",
#   pkg_title = "PGA-cmdr – Population Genetics Shiny Application",
#   pkg_description = "Interactive Shiny application for exploratory and
#     descriptive population genetics analyses from multilocus datasets.",
#   author_first_name = "Naffiou",
#   author_last_name  = "Kadiri",
#   author_email      = "naffiou.kadiri@ird.fr",
#   repo_url          = "https://github.com/KadNaf/PGA1.git"
# )

# 4) Set golem options
golem::set_golem_options()

# 5) Add standard golem files
golem::use_recommended_tests()
golem::use_recommended_deps()

# 6) Compile Rcpp attributes (must be run from pgacmdr/ root)
Rcpp::compileAttributes()

message("01_start.R complete – check DESCRIPTION and NAMESPACE.")
