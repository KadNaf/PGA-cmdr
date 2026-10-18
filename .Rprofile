# Only activate renv in interactive sessions or explicit dev mode.
# Skip during R CMD build / R CMD INSTALL / remotes / pak installs.
if (interactive() || nzchar(Sys.getenv("RENV_ACTIVATE"))) {
  source("renv/activate.R")
}
