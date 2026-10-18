# ── Dockerfile ───────────────────────────────────────────────────────────
# PGA-cmdr - production container
# Build:  docker build -t pgacmdr:latest .
# Run:    docker run --rm -p 3838:3838 pgacmdr:latest

FROM rocker/shiny:latest

LABEL org.opencontainers.image.source="https://github.com/KadNaf/PGA1"
LABEL org.opencontainers.image.licenses="MIT"

# Disable renv and any user .Rprofile for all R commands in this image
ENV R_PROFILE_USER=/dev/null

# ── System libraries ──────────────────────────────────────────────────────
# curl/ssl/xml       : httr, curl, openssl, xml2
# fontconfig/harfbuzz/fribidi/freetype/png/tiff/jpeg : systemfonts, ragg (ggplot2, kableExtra)
# udunits2/gdal/geos/proj : units -> s2 -> sf -> terra/raster -> leaflet
# cmake/ninja-build  : duckdb source compilation fallback on new R versions
# zlib / pandoc      : compression, vignette rendering
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    zlib1g-dev \
    cmake \
    ninja-build \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# ── R package dependencies ────────────────────────────────────────────────
# Use pak for faster installs and better binary resolution
# Mirrors exactly: DESCRIPTION Imports + build-time Suggests
RUN R -e "install.packages('pak', repos = 'https://cloud.r-project.org')" && \
    R -e "pak::pkg_install(c( \
      'shiny', 'bslib', 'shinydashboard', 'shinyWidgets', 'shinyalert', \
      'dplyr', 'tidyr', 'kableExtra', 'zip', \
      'Rcpp', 'tibble', 'waiter', \
      'duckdb', 'DBI', 'jsonlite', 'htmltools', 'magrittr', \
      'knitr', 'rmarkdown', 'testthat', 'digest' \
    ), ask = FALSE)"

# ── Copy package source ───────────────────────────────────────────────────
# renv/, .Rprofile, dev/, *.tar.gz excluded via .dockerignore
WORKDIR /build
COPY . .

# ── Build & install (compiles C++ with OpenMP via src/Makevars) ───────────
RUN R CMD INSTALL --preclean .

# ── Run as a non-root user ─────────────────────────────────────────────────
# rocker/shiny already creates a low-privilege "shiny" user for shiny-server;
# reuse it here so the app never runs as root inside the container. Give it
# ownership of the DuckDB working directory / tmp so writes at runtime don't
# fail once we've dropped root.
RUN mkdir -p /home/shiny/pga-data && chown -R shiny:shiny /home/shiny/pga-data
USER shiny
WORKDIR /home/shiny

# ── Entry point ───────────────────────────────────────────────────────────
EXPOSE 3838
CMD ["R", "-e", "options(shiny.host='0.0.0.0', shiny.port=3838); pgacmdr::run_app()"]