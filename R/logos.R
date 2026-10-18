# logos.R
# Read logo files from the installed package and return base64 data URIs.
# Called once inside app_ui() — no resource path, no external URL dependency.

.read_logo <- function(filename, mime) {
  path <- system.file("app/www", filename, package = "pgacmdr")
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  b64 <- jsonlite::base64_enc(readBin(path, "raw", file.info(path)$size))
  paste0("data:", mime, ";base64,", b64)
}

pga_logo_uris <- function() {
  list(
    ird         = .read_logo("ird_logo.png",         "image/png"),
    cirad       = .read_logo("cirad_logo.png",       "image/png"),
    intertryp   = .read_logo("INTERTRYP_logo.png",   "image/png"),
    ucad        = .read_logo("ucad_logo.png",        "image/png"),
    montpellier = .read_logo("montpellier_logo.png", "image/png")
  )
}