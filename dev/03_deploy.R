# ── 03_deploy.R ──────────────────────────────────────────────────────────
# Deployment helpers for Posit Connect / shinyapps.io / Docker.

# ── Option A: rsconnect (shinyapps.io or Posit Connect) ──────────────────
# rsconnect::deployApp(
#   appDir      = ".",
#   appName     = "pgacmdr",
#   appTitle    = "pgacmdr",
#   forceUpdate = TRUE
# )

# ── Option B: build source package for server installation ────────────────
# devtools::build(".")                   # creates pgacmdr_x.y.z.tar.gz
# install.packages("pgacmdr_x.y.z.tar.gz", repos = NULL, type = "source")

# ── Option C: Docker (see Dockerfile in project root) ─────────────────────
# From the PGA1/ directory:
#   docker build -t pgacmdr:latest .
#   docker run --rm -p 3838:3838 pgacmdr:latest
