# server_LD.R

server_LD <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    
    db_tick <- reactive({ rv$db_tick })
    con_r   <- reactive({ req(rv$con); rv$con })
    
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })
    
    tbl_hf_r <- reactive({
      con <- con_r()
      if (exists("duck_tbl_exists", mode = "function", inherits = TRUE) &&
          exists(".duckdb_get_params", mode = "function", inherits = TRUE) &&
          duck_tbl_exists(con, "params")) {
        p <- .duckdb_get_params(con)
        th <- p$tbl_hf %||% "hf"
        return(as.character(th))
      }
      "hf"
    })
    
    db_ready <- reactive({
      db_tick()
      con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing."),
        shiny::need(DBI::dbExistsTable(con, tbl_hf_r()),   "DuckDB hf table missing.")
      )
      TRUE
    })
    
    base_r <- reactive({
      db_ready()
      
      b <- rv$base_ld %||% rv$base %||% rv$base_r %||% rv$genotype_base
      b <- suppressWarnings(as.integer(b))
      if (length(b) == 1L && is.finite(b) && b > 1L) return(as.integer(b))
      
      con <- con_r()
      if (DBI::dbExistsTable(con, "params") &&
          exists(".duckdb_get_params", mode = "function", inherits = TRUE)) {
        p <- .duckdb_get_params(con)
        b <- suppressWarnings(as.integer(p$base %||% p$base_scalar_full %||% p$base_scalar_preview))
        if (length(b) == 1L && is.finite(b) && b > 1L) return(as.integer(b))
      }
      
      1000L
    })
    
    ld_timing <- reactiveVal(NULL)
    
    ld_data <- reactive({
      db_ready()
      
      con      <- con_r()
      hf_tbl   <- tbl_hf_r()
      meta_tbl <- tbl_meta_r()
      base     <- base_r()
      
      shiny::validate(shiny::need(!is.null(con),      "LD: DuckDB connection is NULL."))
      shiny::validate(shiny::need(!is.null(hf_tbl),   "LD: hf DuckDB table name is NULL."))
      shiny::validate(shiny::need(!is.null(meta_tbl), "LD: meta DuckDB table name is NULL."))
      
      hf_info <- DBI::dbGetQuery(
        con,
        sprintf("PRAGMA table_info(%s)", DBI::dbQuoteIdentifier(con, hf_tbl))
      )
      meta_info <- DBI::dbGetQuery(
        con,
        sprintf("PRAGMA table_info(%s)", DBI::dbQuoteIdentifier(con, meta_tbl))
      )
      
      hf_cols <- hf_info$name
      meta_cols <- meta_info$name
      
      if (all(c("individual", "locus", "g") %in% hf_cols)) {
        hf_ind_col   <- "individual"
        hf_locus_col <- "locus"
        hf_gt_col    <- "g"
      } else if (all(c("indiv_id", "locus_id", "gt") %in% hf_cols)) {
        hf_ind_col   <- "indiv_id"
        hf_locus_col <- "locus_id"
        hf_gt_col    <- "gt"
      } else {
        shiny::validate(
          shiny::need(
            FALSE,
            "LD: hf must contain either (individual,locus,g) or (indiv_id,locus_id,gt)."
          )
        )
      }
      
      if ("individual" %in% meta_cols) {
        meta_ind_col <- "individual"
      } else if ("indiv_id" %in% meta_cols) {
        meta_ind_col <- "indiv_id"
      } else {
        shiny::validate(shiny::need(FALSE, "LD: no individual column found in meta table."))
        }
      
      pop_candidates <- c("Population", "population", "pop", "pop_code")
      pop_col <- pop_candidates[pop_candidates %in% meta_cols][1]
      shiny::validate(shiny::need(!is.na(pop_col), "LD: no population column found in meta table."))
      
      hf_tbl_q   <- as.character(DBI::dbQuoteIdentifier(con, hf_tbl))
      meta_tbl_q <- as.character(DBI::dbQuoteIdentifier(con, meta_tbl))
      hf_ind_q   <- as.character(DBI::dbQuoteIdentifier(con, hf_ind_col))
      hf_locus_q <- as.character(DBI::dbQuoteIdentifier(con, hf_locus_col))
      hf_gt_q    <- as.character(DBI::dbQuoteIdentifier(con, hf_gt_col))
      meta_ind_q <- as.character(DBI::dbQuoteIdentifier(con, meta_ind_col))
      pop_q      <- as.character(DBI::dbQuoteIdentifier(con, pop_col))
      
      # Loci in the ORDER OF THE DATA FILE (first row each locus appears in),
      # NOT alphabetical: the locus order drives the order of the columns and
      # of the locus pairs in the output (B12 X C07, B12 X D12, ...).
      loci_sql <- sprintf("
        SELECT locus FROM (
          SELECT CAST(%s AS VARCHAR) AS locus, MIN(rowid) AS _rank
          FROM %s
          WHERE %s IS NOT NULL
          GROUP BY 1
        ) ORDER BY _rank
      ", hf_locus_q, hf_tbl_q, hf_locus_q)
      
      loci_df <- DBI::dbGetQuery(con, loci_sql)
      loci <- as.character(loci_df$locus)
      shiny::validate(shiny::need(length(loci) > 1, "LD: need at least 2 loci."))
      
      base <- as.integer(base)

      case_exprs <- vapply(loci, function(loc) {
        loc_q   <- as.character(DBI::dbQuoteString(con, loc))
        alias_q <- as.character(DBI::dbQuoteIdentifier(con, loc))
        sprintf("MAX(CASE WHEN locus = %s THEN %s END) AS %s",
                loc_q, hf_gt_q, alias_q)
      }, character(1))

      wide_sql <- sprintf("
        WITH long AS (
          SELECT
            CAST(h.%s AS VARCHAR) AS individual,
            CAST(h.%s AS VARCHAR) AS locus,
            CAST(m.%s AS VARCHAR) AS Population,
            h.%s                  AS gt
          FROM %s h
          LEFT JOIN %s m
            ON CAST(h.%s AS VARCHAR) = CAST(m.%s AS VARCHAR)
          WHERE h.%s IS NOT NULL AND h.%s > 0
        )
        SELECT
          Population,
          individual,
          %s
        FROM long
        GROUP BY Population, individual
        ORDER BY Population, individual
      ",
                          hf_ind_q,
                          hf_locus_q,
                          pop_q,
                          hf_gt_q,
                          hf_tbl_q,
                          meta_tbl_q,
                          hf_ind_q,
                          meta_ind_q,
                          hf_gt_q, hf_gt_q,
                          paste(case_exprs, collapse = ",\n          ")
      )
      
      out <- DBI::dbGetQuery(con, wide_sql)
      
      shiny::validate(shiny::need(nrow(out) > 1, "LD: not enough individuals after reshaping."))
      shiny::validate(shiny::need(ncol(out) > 3, "LD: need at least 2 loci after reshaping."))
      
      fixed_cols <- c("Population", "individual")
      # Keep the locus columns in data order (the order they were built in
      # above) — no sort(): an alphabetical sort here used to scramble them.
      locus_cols <- intersect(loci, setdiff(names(out), fixed_cols))
      out <- out[, c(fixed_cols, locus_cols), drop = FALSE]

      # Populations also in the order of the data file (first appearance in
      # the meta table), not alphabetical.
      pop_order <- DBI::dbGetQuery(con, sprintf("
        SELECT Population FROM (
          SELECT CAST(%s AS VARCHAR) AS Population, MIN(rowid) AS _rank
          FROM %s
          WHERE %s IS NOT NULL
          GROUP BY 1
        ) ORDER BY _rank
      ", pop_q, meta_tbl_q, pop_q))$Population
      out <- out[order(match(out$Population, pop_order), out$individual), , drop = FALSE]
      
      out
    })
    
    loci_names <- reactive({
      df <- ld_data()
      loci <- setdiff(names(df), c("Population", "individual", "Individual"))
      shiny::validate(shiny::need(length(loci) > 1, "Need at least 2 loci to compute LD."))
      loci
    })
    
    # Note: the low-permutations warning that used to pop up on click is now
    # written into the generated parameters file instead (see .run_ld_computation
    # and the parameters writer below) — downloadButton clicks don't have a
    # "before the file starts" interactive step to show a modal in.

    .run_ld_computation <- function() {
      df <- ld_data()
      loci <- loci_names()

      shiny::validate(shiny::need(nrow(df) > 1, "Not enough individuals to compute LD."))
      
      start_time <- Sys.time()
      shinyWidgets::updateProgressBar(session, "LD_progress", value = 10)
      
      geno_mat <- as.matrix(df[, loci, drop = FALSE])
      storage.mode(geno_mat) <- "integer"

      nbperms <- as.integer(input$n_iterations)
      if (is.na(nbperms) || nbperms < 1L) nbperms <- 10000L

      set.seed(1)
      res_cpp <- tryCatch({
        ld_pvalues_cpp(Population = df$Population,
                       geno_mat   = geno_mat,
                       base       = base_r(),
                       nbperms    = nbperms)
      }, error = function(e) {
        showNotification(paste("Error in LD C++ computation:", e$message), type = "error")
        NULL
      })
      if (is.null(res_cpp)) return(NULL)
      
      pv <- data.frame(
        Pair = rownames(res_cpp),
        as.data.frame(res_cpp, check.names = FALSE),
        row.names = NULL,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      shinyWidgets::updateProgressBar(session, "LD_progress", value = 100)
      ld_timing(round(difftime(Sys.time(), start_time, units = "secs"), 1))
      
      pv
    }

    ld_results_store <- reactiveVal(NULL)
    ld_results_reactive <- function() ld_results_store()
    
    
    # ---- Value boxes ----
    output$total_pairs_box <- renderValueBox({
      pv <- ld_results_reactive()
      total <- if (is.null(pv)) 0 else nrow(pv)
      
      valueBox(
        value = total,
        subtitle = HTML("<small>Total locus pairs<br>tested for LD</small>"),
        color = "blue",
        icon = icon("link")
      )
    })
    
    output$significant_pairs_box <- renderValueBox({
      pv <- ld_results_reactive()
      if (is.null(pv) || nrow(pv) == 0) {
        return(valueBox("0 (0%)", HTML("<small>Significant pairs<br>No data</small>"), color = "green", icon = icon("exclamation-triangle")))
      }
      
      min_pvals <- min_pvals_by_pair(pv)
      alpha <- 0.05   # significance threshold used for the summary box
      sig_count <- sum(min_pvals < alpha, na.rm = TRUE)
      pct <- round(100 * sig_count / length(min_pvals), 1)
      
      valueBox(
        value = paste0(sig_count, " (", pct, "%)"),
        subtitle = HTML(paste0("<small>Significant pairs<br>p < ", alpha, "</small>")),
        color = if (sig_count > 0) "yellow" else "green",
        icon = icon("exclamation-triangle")
      )
    })
    
    output$mean_pvalue_box <- renderValueBox({
      pv <- ld_results_reactive()
      if (is.null(pv) || nrow(pv) == 0) {
        return(valueBox("N/A", HTML("<small>Mean p-value<br>No data</small>"), color = "red", icon = icon("calculator")))
      }
      
      all_pvals <- suppressWarnings(as.numeric(unlist(pv[, -1, drop = FALSE])))
      all_pvals <- all_pvals[!is.na(all_pvals) & is.finite(all_pvals)]
      mean_p <- if (length(all_pvals) == 0) NA_real_ else mean(all_pvals)
      
      p_display <- ifelse(is.na(mean_p), "N/A",
                          ifelse(mean_p < 0.001, "< 0.001", format(round(mean_p, 4), nsmall = 4)))
      
      valueBox(
        value = p_display,
        subtitle = HTML("<small>Mean p-value<br>across all tests</small>"),
        color = if (is.na(mean_p)) "red" else if (mean_p < 0.05) "red" else "aqua",
        icon = icon("calculator")
      )
    })
    
    output$analysis_time_ld_box <- renderValueBox({
      req(ld_timing())
      time_sec <- ld_timing()
      time_display <- if (time_sec < 60) paste0(time_sec, " s") else paste0(round(time_sec / 60, 1), " min")
      valueBox(time_display, HTML("<small>Analysis time<br>LD computation</small>"), color = "aqua", icon = icon("clock"))
    })
    
    # -----------------------------#
    # helpers (put inside moduleServer)
    # -----------------------------#
    min_pvals_by_pair <- function(pv) {
      if (is.null(pv) || nrow(pv) == 0) return(numeric(0))
      mat <- as.matrix(pv[, -1, drop = FALSE])
      storage.mode(mat) <- "double"
      v <- apply(mat, 1, min, na.rm = TRUE)
      v[!is.finite(v)] <- 1.0
      v
    }
    
    
    
    # -----------------------------#
    # One button, one click, one action: clicking "Run" IS the download
    # request itself — the LD permutation test runs inside this same
    # content() function before the two files are zipped and streamed back.
    # -----------------------------#
    .write_ld_params <- function(con, pv) {
      np   <- suppressWarnings(as.integer(input$n_iterations))
      loci <- tryCatch(loci_names(), error = function(e) character(0))
      pops <- setdiff(names(pv), c("Pair", "All"))
      hdr <- c(
        "Linkage Disequilibrium - genotypic disequilibrium between all pairs of loci",
        sprintf("Dataset: %s", rv$dataset_filename %||% "default_dataset"),
        sprintf("Loci tested (n = %d): %s", length(loci), paste(loci, collapse = ", ")),
        sprintf("Populations (n = %d): %s", length(pops), paste(pops, collapse = ", ")),
        "Loci, locus pairs and populations are listed in the order of the data file.",
        "",
        "Test: log-likelihood ratio (G) statistic on the genotype x genotype contingency table of each locus pair, within each population.",
        sprintf("Permutation test: %d permutations; genotypes at the second locus are permuted among the individuals of each population.", np),
        "  p-value = (b + 1) / (m + 1), b = number of permuted G >= observed G, m = number of permutations.",
        "  Column 'All': G statistics summed over all populations, compared with the sum of the permuted G statistics of the same permutation round.",
        "Missing data: an individual with a missing genotype at either locus of a pair is left out of that pair's table (pairwise deletion).",
        sprintf("Genotype base: %s", base_r())
      )
      if (!is.na(np) && np < 1000L) {
        hdr <- c(hdr, "",
          "WARNING: fewer than 1000 permutations were requested. The Monte Carlo",
          "p-value formula p = (b + 1) / (m + 1) gives a slight overestimation when m",
          "is small, which may produce unreliable significance calls. A minimum of",
          "1000 permutations is recommended; 10 000 or more for publication-quality",
          "results.")
      }
      writeLines(c(hdr, ""), con = con, useBytes = TRUE)
    }

    # One button, one click, one action: clicking "Run" IS the download
    # request itself — the LD permutation test runs inside this same
    # content() function and the methods block + results are written to ONE
    # plain .txt file (no zip).
    output$ui_ld_out_status <- renderUI({
      tags$p(style = "color:#555;font-size:14px;margin-top:6px;",
        "The results will be saved in ", tags$code(paste0("LD_", Sys.Date(), ".txt")), ".")
    })

    output$run_LD <- downloadHandler(
      filename = function() paste0("LD_", Sys.Date(), ".txt"),
      content  = function(file) {
        pv <- .run_ld_computation()
        ld_results_store(pv)
        req(pv)
        con <- file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        .write_ld_params(con, pv)
        writeLines(sprintf("Linkage disequilibrium p-values (all locus pairs), %s permutations; column All = all populations combined", input$n_iterations), con = con)
        write.table(pv, file = con, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
  })
}