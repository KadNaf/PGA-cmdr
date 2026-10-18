# ── 02_dev.R ─────────────────────────────────────────────────────────────
# Day-to-day development helpers.

# Recompile C++ attributes after editing any .cpp file
Rcpp::compileAttributes()

# Reload the package
devtools::load_all(".")

# Regenerate documentation
devtools::document()

# Run checks
devtools::check()

# Launch the app locally
pgacmdr::run_app()
