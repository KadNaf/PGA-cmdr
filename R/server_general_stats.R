# server_general_stats.R






## =========================================================#
# Server_general_stats ####
## =========================================================#

server_general_stats <- function(id, rv) {
  
  moduleServer(id, function(input, output, session) {
    
    
    # ------------------------------------------------------------------#
    # Parallel controls (server-side defaults)
    # - If you already have UI inputs for threads/seed, replace these.
    # ------------------------------------------------------------------#
    .n_threads <- reactive({
      if (!is.null(input$n_threads) && is.finite(input$n_threads)) {
        return(max(1L, as.integer(input$n_threads)))
        
      }
      # sensible default: use all available cores (or leave one free)
      nc <- NA_integer_
      if (requireNamespace("parallel", quietly = TRUE)) {
        nc <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
      }
      if (!is.finite(nc) || length(nc) != 1L || nc < 1L) nc <- 1L
      max(1L, nc - 1L)
    })
    
    .seed <- reactive({
      # prefer an existing input if you have it
      # global seed set in Import data
      sd <- suppressWarnings(as.numeric(rv$seed))
      if (length(sd) == 1L && is.finite(sd)) return(sd)
      12345
    })
    
    ## =========================================================#
    ## DB truth layer (ONLY sources of truth)
    ## =========================================================#
    
    db_tick <- reactive({ rv$db_tick })
    con_r   <- reactive({ shiny::req(rv$con); rv$con })
    
    tbl_meta_r <- reactive({ rv$tbl_meta %|||% "meta" })
    
    # Prefer params key tbl_hf if you store it; fallback to "hf"
    tbl_hf_r <- reactive({
      con <- con_r()
      if (duck_tbl_exists(con, "params")) {
        p <-  .duckdb_get_params(con)
        th <- p$tbl_hf %|||% "hf"
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
    
    n_pop_db_r <- reactive({
      db_ready()
      con <- con_r()
      tbl <- tbl_meta_r()
      DBI::dbGetQuery(con, sprintf(
        "SELECT COUNT(DISTINCT Population) AS n FROM %s WHERE Population IS NOT NULL",
        sql_ident(con, tbl)
      ))$n[[1]]
    })
    
    params_r <- reactive({
      db_ready()
      .duckdb_get_params(con_r())
    })
    
    base_r <- reactive({
      p <- params_r()
      
      base <- suppressWarnings(as.integer(p$base %|||% p$base_scalar_full %|||% p$base_scalar_preview))
      if (length(base) == 1L && is.finite(base) && base > 1L) return(as.integer(base))
      
      hl <- suppressWarnings(as.integer(p$haplotype_length %|||% p$length_haplotype %|||% p$width_scalar_full))
      shiny::validate(need(length(hl) == 1L && is.finite(hl) && hl > 0L,
                    "params must define base OR haplotype_length/length_haplotype to compute base = 10^L"))
      as.integer(10L ^ hl)
    })
    
    meta_r <- reactive({
      db_ready()
      con <- con_r()
      
      DBI::dbGetQuery(con, sprintf("
    SELECT individual, Population
    FROM %s
    WHERE Population IS NOT NULL
    ORDER BY individual
  ", sql_ident(con, tbl_meta_r())))
    })
    
    hf_mat_r <- reactive({
      db_ready()
      con  <- con_r()
      meta <- meta_r()
      tbl_hf <- tbl_hf_r()
      
      shiny::validate(need(nrow(meta) > 0, "meta table is empty"))
      
      pop_levels <- sort(unique(meta$Population))
      pop_code   <- match(meta$Population, pop_levels)
      
      loci <- DBI::dbGetQuery(con, sprintf(
        "SELECT DISTINCT locus_id FROM %s ORDER BY 1",
        sql_ident(con, tbl_hf)
      ))$locus_id
      shiny::validate(need(length(loci) > 0, "hf table has no loci"))
      
      hf <- DBI::dbGetQuery(con, sprintf("
    SELECT indiv_id AS individual, locus_id, gt AS g
    FROM %s
  ", sql_ident(con, tbl_hf)))
      
      hf$ind_idx   <- match(hf$individual, meta$individual)
      hf$locus_idx <- match(hf$locus_id, loci)
      
      
      N <- nrow(meta); L <- length(loci)
      mat <- matrix(0L, nrow = N, ncol = L + 1L)
      mat[, 1] <- as.integer(pop_code)
      
      ok <- !is.na(hf$ind_idx) & !is.na(hf$locus_idx)
      mat[cbind(hf$ind_idx[ok], hf$locus_idx[ok] + 1L)] <- as.integer(hf$g[ok])

      colnames(mat) <- c("pop", as.character(loci))
      attr(mat, "pop_levels") <- pop_levels

      # -- Reordonner les colonnes selon l'ordre physique DuckDB ----------
      # hf_mat_r() retourne ORDER BY 1 (alphabetique).
      # On reordonne ici une fois pour toutes - tous les caluls en aval
      # heritent du bon ordre sans modification.
      loci_alpha   <- as.character(loci)              # ordre alphabetique actuel
      loci_ordered <- loci_order_r()                  # ordre physique DuckDB

      reorder_idx  <- match(loci_ordered, loci_alpha)
      reorder_idx  <- reorder_idx[!is.na(reorder_idx)]

      # col 1 = pop (inchangee), puis loci dans l'ordre physique
      mat <- mat[, c(1L, reorder_idx + 1L), drop = FALSE]

      # IMPORTANT: matrix column-subsetting with `[` does NOT preserve
      # custom attributes (only dim/dimnames survive) - "pop_levels" was
      # silently dropped by the re-ordering just above, which made the
      # G-test output fall back to generic "Pop1".."PopN" labels. Re-attach it.
      attr(mat, "pop_levels") <- pop_levels

      mat
    })

    # -- Ordre physique des loci depuis DuckDB (MIN(rowid)) ------------------
    # Meme logique que locus_order_cte() dans server_allele_frequencies
    loci_order_r <- reactive({
      db_ready()
      con     <- con_r()
      tbl_hf  <- tbl_hf_r()

      # Detecter le nom de la colonne locus
      info      <- DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)",
                    DBI::dbQuoteIdentifier(con, tbl_hf)))
      locus_col <- if ("locus" %in% info$name) "locus" else "locus_id"
      locus_q   <- as.character(DBI::dbQuoteIdentifier(con, locus_col))
      hf_q      <- as.character(DBI::dbQuoteIdentifier(con, tbl_hf))

      as.character(DBI::dbGetQuery(con, sprintf("
        WITH locus_order AS (
          SELECT CAST(%s AS VARCHAR) AS _lo_marker,
                MIN(rowid)          AS _lo_rank
          FROM %s
          GROUP BY CAST(%s AS VARCHAR)
        )
        SELECT DISTINCT CAST(%s AS VARCHAR) AS Marker, lo._lo_rank
        FROM %s h
        LEFT JOIN locus_order lo
          ON CAST(%s AS VARCHAR) = lo._lo_marker
        ORDER BY lo._lo_rank ASC",
        locus_q, hf_q, locus_q,
        locus_q, hf_q, locus_q
      ))$Marker)
    })
    
    
    ## =========================================================#
    ## Containers
    ## =========================================================#
    
    # Basic stats
    result_stats_reactive <- reactiveVal(NULL)
    result_stats_download <- reactiveVal(NULL)
    result_stats_numeric_reactive <- reactiveVal(NULL)
    result_stats_display_rv <- reactiveVal(NULL)
    
    # Bootstrap FIS Analysis
    fis_boot_results <- reactiveVal(NULL)
    fis_boot_results_pop <- reactiveVal(NULL)
    fis_boot_timing <- reactiveVal(NULL)
    perm_results <- reactiveVal(NULL)
    fis_allele_results <- reactiveVal(NULL)
    
    # For "Overall" stats selection
    result_stats_select_reactive <- reactiveVal(NULL)
    
    ## =========================================================#
    ## Convenience reactives for metadata ####
    ## =========================================================#
    
    hs_by_pop_wide_r <- reactive({
      db_ready()
      con  <- con_r()
      base <- base_r()

      long <- duck_hs_by_pop_locus_long(
        con          = con,
        tbl_hf       = tbl_hf_r(),
        tbl_meta     = tbl_meta_r(),
        base         = base,
        missing_code = 0L
      )

      shiny::validate(need(nrow(long) > 0, "No Hs results available (check hf/meta tables)."))

      wide <- tidyr::pivot_wider(
        long,
        names_from  = Population,
        values_from = Hs
      )

      wide <- as.data.frame(wide, stringsAsFactors = FALSE)

      # -- Reordonner les lignes selon l'ordre physique DuckDB ----------------
      # wide$Locus est en ordre alphabetique (ORDER BY 1 dans duck_hs_by_pop_locus_long)
      # loci_order_r() donne l'ordre physique MIN(rowid)
      loci_ordered <- loci_order_r()
      reorder_idx  <- match(loci_ordered, wide$Locus)
      reorder_idx  <- reorder_idx[!is.na(reorder_idx)]

      if (length(reorder_idx) > 0)
        wide <- wide[reorder_idx, , drop = FALSE]

      wide
    })

    ## =========================================================#
    ## Basic stats cache (DB-first)
    ## =========================================================#
    
    basic_cache <- reactiveVal(list(key = NULL, value = NULL))
    
    .get_basic_stats_cached <- function() {
      db_ready()
      
      mat  <- hf_mat_r()
      base <- base_r()
      k    <- n_pop_db_r()
      ui_hash <- .hash_key(list(
        missing_code = 0L
      ))
      key <- .hash_key(list(
        db_tick = db_tick(),
        tbl_hf  = tbl_hf_r(),
        tbl_meta= tbl_meta_r(),
        params_hash = .hash_key(params_r()),
        user_inputs_hash = ui_hash
      ))
      
      
      cur <- basic_cache()
      if (isTRUE(identical(cur$key, key)) && !is.null(cur$value)) {
        return(cur$value)
      }
      
      res <- .compute_basic_stats(mat = mat, base = base, k = k)
      basic_cache(list(key = key, value = res))
      res
    }
    
    ## =========================================================#
    ## Observer Basic stats ####
    ## =========================================================#
    
    .run_basic_stats_computation <- function() {
      db_ready()
      
      result_stats <- .get_basic_stats_cached()
      
      # keep if any other outputs depend on the full table
      result_stats_reactive(result_stats)
      
      selected_stats <- c(
        "Ho"                 = isTRUE(input$ho_checkbox),
        "Hs"                 = isTRUE(input$hs_checkbox),
        "Ht"                 = isTRUE(input$ht_checkbox),
        # Output order: FIS, FST, FIT (requested convention for all outputs)
        "Fis (W&C)"          = isTRUE(input$fis_wc_checkbox),
        "Fst (W&C)"          = isTRUE(input$fst_wc_checkbox),
        "Fit (W&C)"          = isTRUE(input$fit_wc_checkbox),
        "Fst-max (Meirmans)"  = isTRUE(input$fst_max_checkbox),
        "Fst' (Meirmans)"     = isTRUE(input$fst_prim_checkbox),
        "GST"                = isTRUE(input$GST_checkbox),
        "GST''"              = isTRUE(input$GST_sec_checkbox)
      )
      
      # Only keep columns that exist in result_stats
      keep <- names(selected_stats)[selected_stats]
      
      keep <- intersect(keep, colnames(result_stats))
      
      if (length(keep) > 0) {
        
        result_stats_select <- result_stats[, c("ID", keep), drop = FALSE]
        
        # guarantee Overall last
        if ("Overall" %in% result_stats_select$ID) {
          result_stats_select <- rbind(
            result_stats_select[result_stats_select$ID != "Overall", , drop = FALSE],
            result_stats_select[result_stats_select$ID == "Overall", , drop = FALSE]
          )
        }
        
        result_stats_numeric_reactive(result_stats_select)
        
        result_stats_display <- format_numeric_cols(
          result_stats_select,
          digits = 5,
          exclude = "ID"
        )
        result_stats_download(result_stats_display)
        result_stats_select_reactive(result_stats_select)
        result_stats_display_rv(result_stats_display)
        
        
      } else {
        showNotification("No valid statistics to display", type = "warning")
        result_stats_download(NULL)
        result_stats_numeric_reactive(NULL)
        return(NULL)
      }
      result_stats_select
    }

    
    # -- One button, one click, one action: clicking "Run" IS the download
    #    request itself - computes the selected basic statistics, the
    #    always-available gene-diversity / by-population tables (for ALL
    #    populations automatically), AND the per-allele F-statistics, all
    #    bundled into one zip.
    output$ui_gs_out_status <- renderUI({
      tags$p(style = "color:#555;font-size:14px;margin-top:6px;",
        "The results will be saved in ", tags$code(paste0("general_stats_", Sys.Date(), ".zip")), ".")
    })

    output$run_basic_stats <- downloadHandler(
      filename = function() paste0("general_stats_", Sys.Date(), ".zip"),
      content  = function(file) {
        result_stats_select <- .run_basic_stats_computation()
        req(result_stats_select)

        tmpdir <- tempfile("pga_gs_export_"); dir.create(tmpdir)
        on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

        p1 <- file.path(tmpdir, paste0("basic_statistics_", Sys.Date(), ".txt"))
        write.table(result_stats_select, file = p1, sep = "\t", row.names = FALSE, quote = FALSE)

        p2 <- file.path(tmpdir, paste0("gene_diversity_hs_by_pop_", Sys.Date(), ".txt"))
        hs_w <- hs_by_pop_wide_r()
        obp0 <- pga_overall_by_population(hf_mat_r(), base_r(), missing_code = 0L)
        hs_all <- stats::setNames(rep(NA_real_, ncol(hs_w)), names(hs_w))
        hs_all[[1]] <- NA
        pm <- match(names(hs_w), as.character(obp0$overall[[1]]))
        for (j in seq_along(hs_w)) if (!is.na(pm[j])) hs_all[[j]] <- obp0$overall$Hs[pm[j]]
        hs_w[[1]] <- as.character(hs_w[[1]])
        hs_w <- rbind(hs_w, as.data.frame(as.list(hs_all), check.names = FALSE, stringsAsFactors = FALSE))
        hs_w[[1]][nrow(hs_w)] <- "All loci"
        write.table(hs_w, file = p2, sep = "\t", row.names = FALSE, quote = FALSE)

        # Per-population summary: Ho and Hs averaged over loci weighted by the number of
        # individuals typed at each locus (GENEPOP convention); FIS = Weir & Cockerham.
        obp <- pga_overall_by_population(hf_mat_r(), base_r(), missing_code = 0L)

        # FSTAT-style file: FIS (WC84) per locus x sub-sample, with an "All" column
        # (per locus, over all sub-samples) and an "All" row (per sub-sample, over all
        # loci); bottom-right cell = FIS over all sub-samples and loci.
        p3 <- file.path(tmpdir, paste0("FIS_per_locus_and_subsample_", Sys.Date(), ".txt"))
        con3 <- file(p3, open = "w", encoding = "UTF-8")
        mat3 <- as.matrix(hf_mat_r()); storage.mode(mat3) <- "integer"
        lp3 <- pga_fis_locus_pop(mat3, base = base_r(), pop_names = as.character(obp$overall[[1]]),
                                 n_boot = 0L, n_perm = 0L)
        w3 <- tapply(lp3$FIS, list(factor(lp3$Locus, unique(lp3$Locus)),
                                   factor(lp3$Sub_sample, unique(lp3$Sub_sample))), identity)
        w3 <- matrix(as.numeric(w3), nrow = nlevels(factor(lp3$Locus, unique(lp3$Locus))),
                     dimnames = list(unique(lp3$Locus), unique(lp3$Sub_sample)))
        o3 <- observed_wc84_stats_cpp(mat3, 1L, 0L, as.integer(base_r()))
        pop3 <- tryCatch(as.numeric(wc_fis_by_pop_wc84(mat3, 0L, as.integer(base_r()))), error = function(e) rep(NA_real_, ncol(w3)))
        full3 <- rbind(cbind(w3, All = as.numeric(o3$FIS)),
                       All = c(pop3, as.numeric(o3$FIS_overall_ratio_of_sums)))
        writeLines(c(
          "FIS (Weir & Cockerham 1984) per locus and sub-sample.",
          "Column 'All' = FIS of each locus over all sub-samples; row 'All' = FIS of each sub-sample over all loci; bottom-right cell = FIS over all sub-samples and loci (multilocus ratio of sums, identical to the Overall row of the Local panmixia tables).",
          ""), con = con3, useBytes = TRUE)
        write.table(data.frame(Locus = rownames(full3), full3, check.names = FALSE),
                    file = con3, sep = "\t", row.names = FALSE, quote = FALSE)
        close(con3)

        all_files <- c(p1, p2, p3)

        # Per-allele F-statistics (formerly its own separate Run button) -
        # computed automatically here too, so a single click covers both.
        fres <- .run_allele_fstats_computation()
        if (!is.null(fres)) {
          p4 <- file.path(tmpdir, paste0("Fstats_per_allele_", Sys.Date(), ".txt"))
          write.table(fres, file = p4, sep = "\t", row.names = FALSE, quote = FALSE)
          all_files <- c(all_files, p4)
        }

        zip::zip(zipfile = file, files = basename(all_files), root = tmpdir)
      }
    )
    # =========================================================#
    ### Population-specific stats (DB-native, vectorised) ####
    # =========================================================#
    
    # ---- UI selector choices (no loops, DB query only)
    observe({
      db_ready()
      con <- con_r()
      
      df <- DBI::dbGetQuery(con, sprintf("
    SELECT DISTINCT Population
    FROM %s
    WHERE Population IS NOT NULL
    ORDER BY Population
  ", sql_ident(con, tbl_meta_r())))
      
      shiny::validate(need(nrow(df) > 0, "No populations available in meta table yet."))
      choices <- c("All", df$Population)
      
      updateSelectInput(session, "selected_pop_overall",
                        choices = choices,
                        selected = "All")
    })

    
    # =========================================================#
    ## Download handlers (population section) ####
    # =========================================================#
    
    # 2) Overall by population (Ho/Hs/Fis Nei) - rendered on-screen only,
    #    downloadable via the merged Run+Download button above.

    
    # ==================================== FIS SECTION ANALYSIS ===============================================
    run_fis_by_pop <- function(
    n_perm,
    n_boot,
    conf_level,
    missing_code = 0L
    ) {
      mat  <- hf_mat_r()
      base <- base_r()
      con  <- con_r()
      
      pop_df <- DBI::dbGetQuery(con, sprintf("
    SELECT DISTINCT Population
    FROM %s
    WHERE Population IS NOT NULL
    ORDER BY Population
  ", sql_ident(con, tbl_meta_r())))
      
      pop_codes <- as.character(sort(unique(mat[, 1])))
      pop_names <- as.character(pop_df$Population[seq_along(pop_codes)])
      
      pop_lookup <- stats::setNames(pop_names, pop_codes)
      
      # 2) Observed population-wise FIS - WC84 ratio-of-sums (C++)
      obs_pop <- wc_fis_by_pop_wc84(
        dat     = mat,
        pop_col = 0L,
        base    = base
      )
      obs_codes <- names(obs_pop)
      if (is.null(obs_codes) || !length(obs_codes)) {
        obs_codes <- pop_codes
        names(obs_pop) <- obs_codes
      } else {
        obs_codes <- as.character(obs_codes)
      }
      
      names(obs_pop) <- unname(pop_lookup[obs_codes])
      
      # 3) Bootstrap individuals within populations (C++)
      # This routine draws from R's global random number generator: without a seed its
      # result depended on everything run earlier in the session (not reproducible).
      set.seed(as.integer(.seed()))
      boot_mat <- boot_indiv_wc_fis_by_pop(
        mat     = mat,
        pop_col = 0L,
        NAcode  = as.integer(missing_code),
        B       = as.integer(n_boot),
        base    = base
      )
      # rename bootstrap columns to pop_names
      if (!is.null(colnames(boot_mat))) {
        boot_codes <- as.character(colnames(boot_mat))
        colnames(boot_mat) <- unname(pop_lookup[boot_codes])
      }
      
      alpha <- (1 - conf_level) / 2
      boot_mean <- colMeans(boot_mat, na.rm = TRUE)
      ci_l <- apply(boot_mat, 2, stats::quantile, probs = alpha,     na.rm = TRUE, type = 7)
      ci_u <- apply(boot_mat, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE, type = 7)

      # Fewer than 5 individuals in a given sub-sample (population): its
      # individual bootstrap is not meaningful - report NA for that
      # population specifically, not the whole analysis.
      pop_sizes <- table(as.character(mat[, 1]))
      names(pop_sizes) <- unname(pop_lookup[names(pop_sizes)])
      small_pops <- names(pop_sizes)[pop_sizes < 5L]
      small_pops <- intersect(small_pops, names(boot_mean))
      if (length(small_pops)) {
        boot_mean[small_pops] <- NA_real_
        ci_l[small_pops] <- NA_real_
        ci_u[small_pops] <- NA_real_
      }
      
      # 4) Permutation test (C++) + p-values (two-sided abs)
      perm_res <- NULL
      pvals <- rep(NA_real_, length(pop_names))
      names(pvals) <- pop_names
      
      if (!is.null(n_perm) && n_perm > 0) {
        set.seed(as.integer(.seed()))     # same reason as above (R's global generator)
        perm_res <- batch_permute_wc_fis_by_pop(
          dat            = mat,
          pop_col_1based = 1L,
          base           = as.integer(base),
          B              = as.integer(n_perm)
        )
        
        # rename permutation columns to pop_names
        if (!is.null(colnames(perm_res))) {
          perm_codes <- as.character(colnames(perm_res))
          colnames(perm_res) <- unname(pop_lookup[perm_codes])
        }
        
        common_names <- intersect(names(obs_pop), colnames(perm_res))
        pvals <- rep(NA_real_, length(names(obs_pop)))
        names(pvals) <- names(obs_pop)
        
        pvals[common_names] <- vapply(common_names, function(pn) {
          permj <- perm_res[, pn]
          permj <- permj[is.finite(permj)]
          obsj  <- obs_pop[pn]
          if (!is.finite(obsj) || length(permj) == 0) return(NA_real_)
          ge <- pga_n_abs_ge(permj, obsj)
          (ge + 1) / (length(permj) + 1)
        }, numeric(1))
        
      }
      
      # 5) Overall row (mean across populations)
      obs_overall <- mean(obs_pop, na.rm = TRUE)
      
      boot_overall <- rowMeans(boot_mat, na.rm = TRUE)
      boot_overall <- boot_overall[is.finite(boot_overall)]
      
      if (length(boot_overall) > 0) {
        overall_boot_mean <- mean(boot_overall, na.rm = TRUE)
        overall_ci_l <- as.numeric(stats::quantile(boot_overall, probs = alpha,     na.rm = TRUE, type = 7))
        overall_ci_u <- as.numeric(stats::quantile(boot_overall, probs = 1 - alpha, na.rm = TRUE, type = 7))
      } else {
        overall_boot_mean <- NA_real_
        overall_ci_l <- NA_real_
        overall_ci_u <- NA_real_
      }
      
      overall_p <- NA_real_
      if (!is.null(perm_res)) {
        perm_overall <- rowMeans(perm_res, na.rm = TRUE)
        perm_overall <- perm_overall[is.finite(perm_overall)]
        if (length(perm_overall) > 0 && is.finite(obs_overall)) {
          ge <- pga_n_abs_ge(perm_overall, obs_overall)
          overall_p <- (ge + 1) / (length(perm_overall) + 1)
        }
      }
      
      # 6) Final table
      final_df <- data.frame(
        ID           = pop_names,
        Observed_FIS = as.numeric(obs_pop[pop_names]),
        Boot_Mean    = as.numeric(boot_mean[pop_names]),
        SE           = as.numeric(apply(boot_mat, 2, stats::sd, na.rm = TRUE)[pop_names]),
        CI_L         = as.numeric(ci_l[pop_names]),
        CI_U         = as.numeric(ci_u[pop_names]),
        P_value      = as.numeric(pvals[pop_names]),
        stringsAsFactors = FALSE
      )
      
      overall_row <- data.frame(
        ID           = "Overall",
        Observed_FIS = obs_overall,
        Boot_Mean    = overall_boot_mean,
        SE           = stats::sd(boot_overall, na.rm = TRUE),
        CI_L         = overall_ci_l,
        CI_U         = overall_ci_u,
        P_value      = overall_p,
        stringsAsFactors = FALSE
      )
      
      final_df <- rbind(final_df, overall_row)
      final_df$SE[final_df$ID %in% small_pops] <- NA_real_
      
      list(
        final_table          = final_df,
        permutation_results  = perm_res,
        bootstrap_results    = boot_mat,
        metadata = list(
          id_col     = "Population",
          n_perm     = n_perm,
          n_boot     = n_boot,
          conf_level = conf_level,
          base       = base,
          pop_names  = pop_names
        )
      )
    }
    
    
    run_fis_by_locus <- function(
    n_perm = 1000,
    n_boot = 1000,
    conf_level = 0.95,
    missing_code = 0L
    ) {
      # 1) Source from DuckDB (Design 1)
      mat  <- hf_mat_r()
      base <- base_r()

      # Safety: integer matrix + pop codes in col 1 (R is 1-based)
      mat <- as.matrix(mat)
      storage.mode(mat) <- "integer"
      stopifnot(is.integer(mat))
      stopifnot(ncol(mat) >= 2L)
      stopifnot(all(mat[, 1] > 0, na.rm = TRUE))
      stopifnot(!is.na(base), base > 0L)
      
      # Locus names: all columns except population code column
      locus_names <- colnames(mat)[-1]
      if (is.null(locus_names) || !length(locus_names)) {
        locus_names <- paste0("L", seq_len(ncol(mat) - 1L))
      }
      
      # 2) Observed FIS (C++)
      observed_fis <- fis_wc_cpp(mat, base = as.integer(base))$FIS
      # 2b) Correct overall FIS - WC84 ratio-of-sums: sum(B) / sum(B+C) across loci
      # This is NOT the mean of per-locus FIS values.
      fis_overall_obs_stats <- observed_wc84_stats_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      )
      fis_obs_overall <- as.numeric(fis_overall_obs_stats$FIS_overall_ratio_of_sums)
      # 3) Bootstrap individuals within populations (C++)
      boot_mat <- boot_indiv_wc_fis(
        mat     = mat,
        pop_col = 0L,
        NAcode  = as.integer(missing_code),
        B       = as.integer(n_boot),
        base    = as.integer(base)
      )
      # 4) Bootstrap summaries (individual bootstrap - always computed)
      summary_list <- summarize_fis_results(
        boot = boot_mat,
        conf = conf_level
      )
      # 3b/4b) Bootstrap POPULATIONS (pop-block) - captures uncertainty from
      # the sampling of populations, not just individuals within populations.
      # This is the same resampling scheme used for FST.
      n_subs_avail_fis <- length(unique(mat[, 1]))
      if (n_subs_avail_fis < 5L) {
        # Fewer than 5 sub-samples (populations): bootstrap over sub-samples
        # is not meaningful - report NA instead of computing.
        na_loc <- rep(NA_real_, length(locus_names))
        summary_pop_list <- list(
          mean = na_loc, ci_lower = na_loc, ci_upper = na_loc,
          overall_mean = NA_real_, overall_ci_lower = NA_real_, overall_ci_upper = NA_real_
        )
      } else {
        boot_pop_mat <- boot_popblock_wc_fis(
          mat     = mat,
          pop_col = 0L,
          NAcode  = as.integer(missing_code),
          B       = as.integer(n_boot),
          base    = as.integer(base)
        )
        summary_pop_list <- summarize_fis_results(
          boot = boot_pop_mat,
          conf = conf_level
        )
      }
      # 5) Build final dataframe (Observed + CI from individual bootstrap, matching Genetix)
      final_df <- create_results_dataframe(
        obs         = observed_fis,
        sum         = summary_list,
        locus_names = locus_names
      )

      # ---- Canonicalise identifier column to ID
      id_src <- NULL
      if ("Locus" %in% names(final_df)) {
        id_src <- "Locus"
      } else if ("ID" %in% names(final_df)) {
        id_src <- "ID"
      } else if (length(locus_names) == nrow(final_df) - 1L) {
        final_df <- tibble::rownames_to_column(final_df, var = "Locus")
        id_src <- "Locus"
      } else {
        stop(
          "run_fis_by_locus(): final_df has no identifier column. Columns are: ",
          paste(names(final_df), collapse = ", ")
        )
      }
      final_df <- dplyr::rename(final_df, ID = dplyr::all_of(id_src))
      # Add Boot_Mean from individual bootstrap.
      final_df$Boot_Mean <- c(
        as.numeric(summary_list$mean),
        as.numeric(summary_list$overall_mean)
      )
      # Override the Overall row with the correct ratio-of-sums observed and CI.
      overall_row_idx <- which(final_df$ID == "Overall")
      if (length(overall_row_idx) == 1L) {
        final_df$Observed_FIS[overall_row_idx] <- fis_obs_overall
        final_df$CI_L[overall_row_idx]         <- summary_list$overall_ci_lower
        final_df$CI_U[overall_row_idx]         <- summary_list$overall_ci_upper
        final_df$Boot_Mean[overall_row_idx]    <- summary_list$overall_mean
      }
      # 6) Permutation test p-values (two-sided abs)
      perm_res <- NULL
      
      if (!is.null(n_perm) && n_perm > 0) {
        
        perm_res <- batch_permute_wc_fis(
          dat            = mat,
          pop_col_1based = 1L,
          base           = as.integer(base),
          B              = as.integer(n_perm)
        )
        
        perm_res_loci <- perm_res[, seq_len(ncol(perm_res) - 1L), drop = FALSE]

        pvals_loci <- vapply(seq_along(observed_fis), function(j) {
          permj <- perm_res_loci[, j]
          permj <- permj[is.finite(permj)]
          obsj  <- observed_fis[j]
          if (!is.finite(obsj) || length(permj) == 0) return(NA_real_)
          ge <- pga_n_abs_ge(permj, obsj)
          (ge + 1) / (length(permj) + 1)
        }, numeric(1))

        # Last column of perm_res is the ratio-of-sums overall FIS per replicate
        perm_overall <- perm_res[, ncol(perm_res)]
        perm_overall <- perm_overall[is.finite(perm_overall)]
        
        p_overall <- if (is.finite(fis_obs_overall) && length(perm_overall) > 0) {
          ge <- pga_n_abs_ge(perm_overall, fis_obs_overall)
          (ge + 1) / (length(perm_overall) + 1)
        } else {
          NA_real_
        }
        
        final_df$P_value <- c(pvals_loci, p_overall)
        
        
      } else {
        final_df$P_value <- rep(NA_real_, nrow(final_df))
      }
      
      list(
        final_table             = final_df,
        permutation_results     = perm_res,
        bootstrap_results       = boot_mat,        # individual bootstrap (within pops)
        bootstrap_pop_results   = boot_pop_mat,    # pop-block bootstrap
        ci_indiv = list(          # CIs from individual bootstrap
          mean  = summary_list$mean,
          ci_lo = summary_list$ci_lower,
          ci_hi = summary_list$ci_upper,
          overall_mean  = summary_list$overall_mean,
          overall_ci_lo = summary_list$overall_ci_lower,
          overall_ci_hi = summary_list$overall_ci_upper
        ),
        ci_pop = list(            # CIs from pop-block bootstrap
          mean  = summary_pop_list$mean,
          ci_lo = summary_pop_list$ci_lower,
          ci_hi = summary_pop_list$ci_upper,
          overall_mean  = summary_pop_list$overall_mean,
          overall_ci_lo = summary_pop_list$overall_ci_lower,
          overall_ci_hi = summary_pop_list$overall_ci_upper
        ),
        metadata = list(
          id_col      = "Locus",
          n_perm      = n_perm,
          n_boot      = n_boot,
          conf_level  = conf_level,
          base        = as.integer(base),
          locus_names = locus_names
        )
      )
    }
    
    
    .run_fis_computation <- function() {
      db_ready()
      
      if (input$n_perm < 10 || input$n_boot < 10) {
        showNotification(
          "Number of permutations and bootstrap replicates should be at least 10",
          type = "warning"
        )
      }
      
      results <- tryCatch({
        start_time <- Sys.time()
        shinyWidgets::updateProgressBar(session, "fis_progress", value = 5)
        
        # Both analysis levels are always computed - no selector needed.
        results_locus <- run_fis_by_locus(
          n_perm       = input$n_perm,
          n_boot       = input$n_boot,
          conf_level   = input$conf_level,
          missing_code = 0L
        )
        shinyWidgets::updateProgressBar(session, "fis_progress", value = 50)
        results_pop <- run_fis_by_pop(
          n_perm       = input$n_perm,
          n_boot       = input$n_boot,
          conf_level   = input$conf_level,
          missing_code = 0L
        )
        
        # Bootstrap over LOCI (requested for the output file): global multilocus FIS, and
        # FIS of each population (see pga_boot_over_loci / pga_fis_pop_loci_boot).
        mat_fis  <- hf_mat_r(); base_fis <- as.integer(base_r())
        comp_fis <- wc84_locus_components_cpp(dat = mat_fis, pop_col_1based = 1L,
                                              missing_code = 0L, base = base_fis)
        loci_boot <- pga_boot_over_loci(comp_fis$A, comp_fis$B, comp_fis$C, comp_fis$HS, comp_fis$HT,
                                        n_boot = input$n_boot, conf_level = input$conf_level, seed = .seed())
        pop_loci_boot <- pga_fis_pop_loci_boot(mat_fis, base_fis, input$n_boot, input$conf_level, seed = .seed())

        shinyWidgets::updateProgressBar(session, "fis_progress", value = 100)

        fis_boot_timing(round(difftime(Sys.time(), start_time, units = "secs"), 1))
        # The on-screen value boxes stay driven by the by-locus view.
        fis_boot_results(results_locus)
        fis_boot_results_pop(results_pop)
        perm_results(results_locus$permutation_results)

        # Per-allele F-stats (FIS, FST, FIT via WC84 components)
        allele_mat  <- as.matrix(hf_mat_r())
        storage.mode(allele_mat) <- "integer"
        allele_base <- as.integer(base_r())
        fis_allele_results(wc84_per_allele_fstats_cpp(
          dat          = allele_mat,
          pop_col      = 0L,
          base         = allele_base,
          missing_code = 0L
        ))

        showNotification("Computations completed.", type = "message")
        list(by_locus = results_locus, by_pop = results_pop,
             loci_boot = loci_boot, pop_loci_boot = pop_loci_boot)

      }, error = function(e) {
        fis_boot_results(NULL)
        fis_boot_results_pop(NULL)
        perm_results(NULL)
        fis_allele_results(NULL)
        showNotification(paste("Error in bootstrap analysis:", e$message), type = "error")
        NULL
      })

      results
    }

    
    ## FIS value boxes ----
    
    
    
    
    
    
    
    # -- Methods/parameters block written at the TOP of the single FIS output file.
    .write_fis_params <- function(con, res) {
      nb <- input$n_boot
      hdr <- c(
        "Local Panmixia - FIS (Weir & Cockerham 1984)",
        sprintf("Dataset: %s", rv$dataset_filename %||% "default_dataset"),
        "",
        sprintf("P-values: permutation test, %s permutations; alleles randomised among individuals within each sub-sample.", input$n_perm),
        "  P_value (two-sided) = (b + 1) / (m + 1), b = number of permuted |FIS| >= observed |FIS|, m = number of permutations.",
        "  Table 1 also gives the two one-sided p-values of FSTAT: P_deficit (heterozygote deficit, permuted FIS >= observed) and P_excess (heterozygote excess, permuted FIS <= observed).",
        "",
        sprintf("Confidence intervals: %s level, percentile bootstrap, %s replicates for each type of bootstrap:", input$conf_level, nb),
        "  _indiv = bootstrap over INDIVIDUALS (individuals resampled within each sub-sample). SE = standard deviation of the bootstrap replicates.",
        "  _subs  = bootstrap over SUB-SAMPLES (sub-samples resampled as blocks); Table 3 only; NA if fewer than 5 sub-samples.",
        "  _loci  = bootstrap over LOCI (loci resampled with replacement); Tables 2, 3 and 4.",
        "  NA rule: a sub-sample with fewer than 5 individuals gets NA for its _indiv interval and SE.",
        "  Overall rows (Tables 2, 3, 4): multilocus ratio-of-sums FIS over all sub-samples and loci (the same value, CI and p-value in the three tables).",
        ""
      )
      writeLines(hdr, con = con, useBytes = TRUE)
    }

    # By-locus table: bootstrap over individuals and over sub-samples for every locus,
    # plus bootstrap over loci for the Overall row.
    .fis_locus_export <- function(res) {
      r  <- res$by_locus
      ft <- r$final_table
      n  <- nrow(ft)
      pick <- function(x, ov) {
        v <- c(as.numeric(x), as.numeric(ov))
        if (length(v) == n) v else rep(NA_real_, n)
      }
      lb <- res$loci_boot$table
      lb <- lb[lb$Statistic == "FIS", , drop = FALSE]
      cl <- ch <- rep(NA_real_, n)
      if (nrow(lb) == 1L) { cl[n] <- lb$CI_L; ch[n] <- lb$CI_U }   # Overall is the last row
      data.frame(
        ID           = ft$ID,
        Observed_FIS = ft$Observed_FIS,
        SE_indiv     = {
          # bootstrap_results has one column per locus and a last "Overall" column
          se <- apply(r$bootstrap_results, 2, stats::sd, na.rm = TRUE)
          if (length(se) == n) as.numeric(se) else rep(NA_real_, n)
        },
        CI_L_indiv   = ft$CI_L,
        CI_U_indiv   = ft$CI_U,
        CI_L_subs    = pick(r$ci_pop$ci_lo, r$ci_pop$overall_ci_lo),
        CI_U_subs    = pick(r$ci_pop$ci_hi, r$ci_pop$overall_ci_hi),
        CI_L_loci    = cl,
        CI_U_loci    = ch,
        P_value      = ft$P_value,
        stringsAsFactors = FALSE
      )
    }

    # By-sub-sample table: bootstrap over individuals and over loci. The "Overall"
    # row is NOT the mean of the sub-sample values: it is the multilocus
    # ratio-of-sums FIS over all sub-samples and loci, i.e. the Overall row of the
    # by-locus table (same observed value, same CI, same p-value).
    .fis_pop_export <- function(res) {
      ft <- res$by_pop$final_table
      pl <- res$pop_loci_boot$table
      m  <- match(ft$ID, pl$Population)
      se <- if ("SE" %in% names(ft)) ft$SE else rep(NA_real_, nrow(ft))
      out <- data.frame(
        ID           = ft$ID,
        Observed_FIS = ft$Observed_FIS,
        SE_indiv     = se,
        CI_L_indiv   = ft$CI_L,
        CI_U_indiv   = ft$CI_U,
        CI_L_loci    = pl$CI_L_loci[m],
        CI_U_loci    = pl$CI_U_loci[m],
        P_value      = ft$P_value,
        stringsAsFactors = FALSE
      )
      lo <- .fis_locus_export(res); ov <- lo[nrow(lo), ]
      n <- nrow(out)
      out$ID[n]           <- "Overall"
      out$Observed_FIS[n] <- ov$Observed_FIS
      out$SE_indiv[n]     <- ov$SE_indiv
      out$CI_L_indiv[n]   <- ov$CI_L_indiv
      out$CI_U_indiv[n]   <- ov$CI_U_indiv
      out$CI_L_loci[n]    <- ov$CI_L_loci
      out$CI_U_loci[n]    <- ov$CI_U_loci
      out$P_value[n]      <- ov$P_value
      out
    }

    # -- One button, one click, one action: clicking "Run" IS the download
    #    request itself - both analysis levels run inside this same
    #    content() function and everything (methods block + By Locus +
    #    By Population) is written to ONE plain .txt file (no zip, no figure).
    output$ui_fis_out_status <- renderUI({
      tags$p(style = "color:#555;font-size:14px;margin-top:6px;",
        "The results will be saved in ", tags$code(paste0("local_panmixia_FIS_", Sys.Date(), ".txt")), ".")
    })

    output$Run_FIS_Analysis <- downloadHandler(
      filename = function() paste0("local_panmixia_FIS_", Sys.Date(), ".txt"),
      content  = function(file) {
        res <- .run_fis_computation()
        req(res)
        # Table 1: FIS per locus x sub-sample (CI, p-values)
        mat_lp <- as.matrix(hf_mat_r()); storage.mode(mat_lp) <- "integer"
        pn <- tryCatch({
          pd <- DBI::dbGetQuery(con_r(), sprintf(
            "SELECT DISTINCT Population FROM %s WHERE Population IS NOT NULL ORDER BY Population",
            sql_ident(con_r(), tbl_meta_r())))
          as.character(pd$Population)[seq_along(sort(unique(mat_lp[, 1])))]
        }, error = function(e) NULL)
        lp <- pga_fis_locus_pop(mat_lp, base = base_r(), pop_names = pn, n_boot = input$n_boot,
                                n_perm = input$n_perm, conf_level = input$conf_level, seed = .seed())
        nb <- input$n_boot; npm <- input$n_perm
        con <- file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        .write_fis_params(con, res)
        lo <- .fis_locus_export(res)
        wt <- function(x) pga_write_table(x, file = con)
        writeLines(sprintf("Table 1: FIS per locus and sub-sample - CI from %s bootstrap over individuals; p-values from %s permutations", nb, npm), con = con)
        wt(lp)
        writeLines("", con = con)
        writeLines(sprintf("Table 2: FIS per sub-sample over all loci - CI from %s bootstrap over individuals and over loci; two-sided p-values from %s permutations", nb, npm), con = con)
        wt(.fis_pop_export(res))
        writeLines("", con = con)
        writeLines(sprintf("Table 3: FIS per locus over all sub-samples - CI from %s bootstrap over individuals, over sub-samples and (Overall row) over loci; two-sided p-values from %s permutations", nb, npm), con = con)
        wt(lo[-nrow(lo), , drop = FALSE])
        writeLines("", con = con)
        writeLines(sprintf("Table 4: FIS over all sub-samples and loci - CI from %s bootstrap over individuals, over sub-samples and over loci; two-sided p-value from %s permutations", nb, npm), con = con)
        ov <- lo[nrow(lo), , drop = FALSE]; ov$ID <- "All sub-samples and loci"
        wt(ov)
      }
    )

    ## Compute per-allele F-statistics (independent of bootstrap analyses) ----
    .run_allele_fstats_computation <- function() {
      shiny::req(hf_mat_r(), base_r())
      res <- tryCatch({
        allele_mat  <- as.matrix(hf_mat_r())
        storage.mode(allele_mat) <- "integer"
        allele_base <- as.integer(base_r())
        out <- wc84_per_allele_fstats_cpp(
          dat          = allele_mat,
          pop_col      = 0L,
          base         = allele_base,
          missing_code = 0L
        )
        # Sort by Locus (in the ORIGINAL data column order - not
        # alphabetically) then Allele (numeric) - the C++ routine doesn't
        # guarantee output order.
        loc_col <- intersect(c("Locus", "locus", "Marker"), names(out))[1]
        all_col <- intersect(c("Allele", "allele"), names(out))[1]
        if (!is.na(loc_col) && !is.na(all_col)) {
          data_locus_order <- colnames(allele_mat)[-1]
          loc_rank <- match(as.character(out[[loc_col]]), data_locus_order)
          ord <- order(loc_rank, suppressWarnings(as.numeric(out[[all_col]])))
          out <- out[ord, , drop = FALSE]
        }
        fis_allele_results(out)
        showNotification("Per-allele F-statistics computed successfully!", type = "message")
        out
      }, error = function(e) {
        fis_allele_results(NULL)
        showNotification(paste("Error computing per-allele F-statistics:", e$message), type = "error")
        NULL
      })
      res
    }

    # ==================================== FIT SECTION ANALYSIS ===============================================
    ## FIT Analysis reactives ----
    fit_boot_results <- reactiveVal(NULL)
    fit_boot_timing  <- reactiveVal(NULL)
    fit_perm_results <- reactiveVal(NULL)
    ## FIT boot and perm analysis ----
    run_fit_analysis <- function(
    n_perm = 1000,
    n_boot = 1000,
    conf_level = 0.95,
    missing_code = 0L,
    cpp_file = "src/fit_permute_bootstrap_wc.cpp",
    cpp_verbose = FALSE
    ) {
      
      
      # DB-first source
      db_ready()
      mat  <- hf_mat_r()
      base <- base_r()
      
      mat <- as.matrix(mat)
      storage.mode(mat) <- "integer"
      shiny::validate(
        shiny::need(is.integer(mat), "hf_mat_r() must return an integer matrix"),
        shiny::need(ncol(mat) >= 2L, "hf_mat_r() must be pop   >= 1 locus"),
        shiny::need(all(mat[, 1] > 0, na.rm = TRUE), "Population codes must be positive integers (1..K)"),
        shiny::need(isTRUE(is.finite(base)) && base > 1L, "base_r() returned invalid base")
      )
      
      locus_names <- colnames(mat)[-1L]
      if (is.null(locus_names) || length(locus_names) == 0L) {
        locus_names <- paste0("L", seq_len(ncol(mat) - 1L))
      }
      
      # observed_wc84_stats_cpp: WC84 per-locus FST/FIT/HI/HS for the FIT bootstrap pipeline.
      # (observed_wc84_stats_cpp from the OpenMP file is the full-stat version used by the FST section.)
      obs_stats <- observed_wc84_stats_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      )

      observed_fit <- obs_stats$FIT
      names(observed_fit) <- locus_names
      
      # BS/Indiv: individuals are resampled WITHIN each sub-sample (sub-sample sizes
      # and the set of sub-samples stay fixed). FIT is unbiased under this scheme,
      # so the plain percentile CI is used (no shift needed, but the shift is
      # applied for consistency: it is ~0 for FIT).
      boot_res <- boot_indiv_wc84_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base),
        B              = as.integer(n_boot),
        seed           = as.numeric(.seed()),
        n_threads      = 0L
      )
      fit_boot_mat <- boot_res$FIT_boot
      colnames(fit_boot_mat) <- locus_names
      ci_mat <- vapply(seq_along(locus_names), function(j)
        pga_ci_bias_shift(observed_fit[j], fit_boot_mat[, j], conf_level), numeric(2))
      ci_lower <- as.numeric(ci_mat[1, ])
      ci_upper <- as.numeric(ci_mat[2, ])
      CI_FIT <- rbind(lo = ci_lower, hi = ci_upper)
      boot_mean <- as.numeric(colMeans(fit_boot_mat, na.rm = TRUE))
      boot_se   <- apply(fit_boot_mat, 2, stats::sd, na.rm = TRUE)
      perm_res <- NULL
      perm_fit_mat <- NULL
      
      pvals_fit <- rep(NA_real_, length(locus_names))
      names(pvals_fit) <- locus_names
      
      if (!is.null(n_perm) && n_perm > 0) {
        
        # FSTAT's "Randomising alleles overall samples" - global allele shuffle
        perm_res <- batch_permute_fit_global(
          dat            = mat,
          pop_col_1based = 1L,
          missing_code   = as.integer(missing_code),
          base           = as.integer(base),
          B              = as.integer(n_perm)
        )
        
        perm_fit_mat <- perm_res$FIT_perm
        colnames(perm_fit_mat) <- locus_names
        
        pvals_fit <- perm_res$p_FIT
        names(pvals_fit) <- locus_names
      }
      
      res_loci <- data.frame(
        ID          = locus_names,
        Observed_FIT = as.numeric(observed_fit),
        Boot_Mean    = boot_mean,
        SE           = as.numeric(boot_se),
        CI_L         = ci_lower,
        CI_U         = ci_upper,
        P_value      = as.numeric(pvals_fit),
        stringsAsFactors = FALSE
      )
      
      # ---- Overall (FIT) ----
      alpha <- (1 - conf_level) / 2
      
      # WC84 ratio-of-sums overall FIT - NOT the mean of per-locus FIT values
      overall_obs <- as.numeric(obs_stats$FIT_overall_ratio_of_sums)
      
      # BS/Indiv overall distribution (ratio-of-sums recomputed on every replicate)
      boot_overall <- boot_res$FIT_overall_boot
      boot_overall <- boot_overall[is.finite(boot_overall)]
      overall_boot_mean <- mean(boot_overall)
      overall_se <- stats::sd(boot_overall)
      if (length(boot_overall) < 10) {
        overall_ci_l <- NA_real_; overall_ci_u <- NA_real_
      } else {
        ci_o <- pga_ci_bias_shift(overall_obs, boot_overall, conf_level)
        overall_ci_l <- ci_o[["lo"]]; overall_ci_u <- ci_o[["hi"]]
      }
      
      # overall p-value from permutation (two-sided abs), if available
      overall_p <- NA_real_
      if (!is.null(perm_res)) {
        # p-value already computed in C++ using ratio-of-sums overall FIT
        overall_p <- as.numeric(perm_res$p_FIT_overall)
      }
      
      overall_row <- data.frame(
        ID          = "Overall",
        Observed_FIT = overall_obs,
        Boot_Mean    = overall_boot_mean,
        SE           = overall_se,
        CI_L         = overall_ci_l,
        CI_U         = overall_ci_u,
        P_value      = overall_p,
        stringsAsFactors = FALSE
      )
      
      
      final_df <- rbind(res_loci, overall_row)
      
      list(
        final_table = final_df,
        observed_fit = observed_fit,
        permutation_results = perm_fit_mat,
        bootstrap_results = list(
          FIT_boot = fit_boot_mat,
          CI_FIT   = CI_FIT
        ),
        metadata = list(
          n_loci = length(locus_names),
          n_permutations = n_perm,
          n_bootstrap = n_boot,
          conf_level = conf_level,
          base = base,
          locus_names = locus_names,
          pval_method_FIT = if (!is.null(perm_res)) perm_res$pval_method_FIT else NA_character_
        )
      )
    }
    
    ## FIT observeEvent bootstrap and permutation (button) ====
    
    .run_fit_computation <- function() {
      db_ready()
      
      if (input$n_perm_fit < 10 || input$n_boot_fit < 10) {
        showNotification(
          "Number of permutations and bootstrap replicates should be at least 10",
          type = "warning"
        )
      }
      
      results <- tryCatch({
        
        start_time <- Sys.time()
        shinyWidgets::updateProgressBar(session, "fit_progress", value = 5)
        
        shinyWidgets::updateProgressBar(session, "fit_progress", value = 12)
        
        # Run analysis
        results <- run_fit_analysis(
          n_perm     = input$n_perm_fit,
          n_boot     = input$n_boot_fit,
          conf_level = input$conf_level_fit,
          missing_code   = 0L,
          cpp_file       = "src/fit_permute_bootstrap_wc.cpp",
          cpp_verbose    = FALSE
        )
        
        shinyWidgets::updateProgressBar(session, "fit_progress", value = 100)
        
        fit_boot_timing(round(difftime(Sys.time(), start_time, units = "secs"), 1))
        fit_boot_results(results)
        
        # If you use perm_results() elsewhere for plotting, keep this:
        fit_perm_results(results$permutation_results)
        
        showNotification("Computations completed.", type = "message")
        results
        
      }, error = function(e) {
        showNotification(paste("Error in FIT analysis:", e$message), type = "error")
        fit_boot_results(NULL)
        fit_perm_results(NULL)
        NULL
      })

      results
    }

    
    ## ===== FIT value boxes =====
    
    
    
    
    
    ## FIT export (single plain .txt: methods block + results; no figure) ====
    .write_fit_params <- function(con, res) {
      hdr <- c(
        "Global Panmixia - FIT (Weir & Cockerham 1984)",
        sprintf("Dataset: %s", rv$dataset_filename %||% "default_dataset"),
        "",
        sprintf("P-values: two-sided permutation test on |FIT|, %s permutations; alleles randomised over all samples (FSTAT 'randomising alleles overall samples').", input$n_perm_fit),
        "  p = (b + 1) / (m + 1), b = number of permuted |FIT| >= observed |FIT|, m = number of permutations.",
        "",
        sprintf("Confidence intervals: %s level, percentile bootstrap, %s replicates for each type of bootstrap:", input$conf_level_fit, input$n_boot_fit),
        "  Locus rows and 'Overall (BS/Indiv)': bootstrap over INDIVIDUALS (individuals resampled with replacement within each sub-sample; the sub-samples are kept). SE = standard deviation of the bootstrap replicates.",
        "  'Overall (BS/Loci)': bootstrap over LOCI (loci resampled with replacement).",
        ""
      )
      writeLines(hdr, con = con, useBytes = TRUE)
    }

    # -- One button, one click, one action: clicking "Run" IS the download
    #    request itself - the FIT bootstrap/permutation runs inside this
    #    same content() function and the methods block + results are written
    #    to ONE plain .txt file (no zip, no figure).
    output$ui_fit_out_status <- renderUI({
      tags$p(style = "color:#555;font-size:14px;margin-top:6px;",
        "The results will be saved in ", tags$code(paste0("global_panmixia_FIT_", Sys.Date(), ".txt")), ".")
    })

    output$Run_FIT_Analysis <- downloadHandler(
      filename = function() paste0("global_panmixia_FIT_", Sys.Date(), ".txt"),
      content  = function(file) {
        res <- .run_fit_computation()
        req(res)
        ft <- res$final_table
        n  <- nrow(ft)                      # loci..., then the Overall row (last)

        # Overall FIT by bootstrap over LOCI (loci resampled with replacement)
        mat_fit <- hf_mat_r()
        comp <- wc84_locus_components_cpp(dat = mat_fit, pop_col_1based = 1L, missing_code = 0L,
                                          base = as.integer(base_r()))
        bl <- pga_boot_over_loci(comp$A, comp$B, comp$C, comp$HS, comp$HT,
                                 n_boot = input$n_boot_fit, conf_level = input$conf_level_fit,
                                 seed = .seed())$table
        bl <- bl[bl$Statistic == "FIT", , drop = FALSE]

        out <- data.frame(
          ID           = ft$ID,
          Observed_FIT = ft$Observed_FIT,
          SE           = ft$SE,
          P_value      = ft$P_value,
          CI_L         = ft$CI_L,
          CI_U         = ft$CI_U,
          stringsAsFactors = FALSE
        )
        out$ID[n] <- "Overall (BS/Indiv)"
        out <- rbind(out, data.frame(ID = "Overall (BS/Loci)", Observed_FIT = ft$Observed_FIT[n],
                                     SE = if ("SE" %in% names(bl)) bl$SE else NA_real_,
                                     P_value = ft$P_value[n], CI_L = bl$CI_L, CI_U = bl$CI_U,
                                     stringsAsFactors = FALSE))

        con <- file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        .write_fit_params(con, res)
        writeLines(sprintf("FIT with confidence intervals from %s bootstrap over individuals (locus rows, Overall BS/Indiv) and over loci (Overall BS/Loci); two-sided p-values from %s permutations",
                           input$n_boot_fit, input$n_perm_fit), con = con)
        pga_write_table(out, file = con)
      }
    )

        # ==================================== FST SECTION ANALYSIS ===============================================
    
    # hard guard: fail early if C++ funcs are missing
    stopifnot(
      exists("observed_wc84_stats_cpp",   mode = "function"),  # OpenMP full-stat version (FST section)
      exists("batch_permute_wc84_fst_parallel", mode = "function"),
      exists("boot_popblock_wc84_parallel",     mode = "function")
    )
    
    
    ## FST Analysis reactives ----
    fst_boot_results <- reactiveVal(NULL)
    fst_boot_timing  <- reactiveVal(NULL)
    fst_parallel_meta <- reactiveVal(NULL)
    # Genetic Diversities used to share fst_boot_results/fst_boot_timing with
    # Subdivision (both call the same run_bootstrap_fst_analysis() helper,
    # which computes FST/HS/HT together) - this silently showed one module's
    # on-screen summary boxes with the OTHER module's last-run results
    # whenever they were run with different parameters. Each module's own
    # downloaded file was never affected (each already computed and wrote
    # its own fresh result directly) - only the value boxes were shared.
    # Giving Diversities its own reactiveVal fixes this without changing
    # what either module computes.
    div_boot_results <- reactiveVal(NULL)
    div_boot_timing  <- reactiveVal(NULL)
    
    # diversity_extras = TRUE  -> also compute the HS-by-individuals bootstrap and the
    #   per-population HS bootstrap (needed ONLY by the Genetic Diversities module).
    # diversity_extras = FALSE -> skipped (Subdivision does not use them; they used to
    #   be computed anyway and made the FST run noticeably slower than it needs to be).
    run_bootstrap_fst_analysis <- function(n_perm, n_boot, conf_level, missing_code = 0L,
                                            progress_id = NULL, diversity_extras = TRUE) {
      db_ready()

      # Real intermediate progress: each call below corresponds to one of the
      # (sequential, often the slowest) stages of the analysis, so updating
      # here reflects actual work done rather than a single 15% -> 100% jump.
      .bump_progress <- function(value, title = NULL) {
        if (!is.null(progress_id)) {
          shinyWidgets::updateProgressBar(session, progress_id, value = value,
                                           title = title %||% "Overall Progress")
        }
      }

      .bump_progress(5, "Loading genotype matrix...")
      mat  <- hf_mat_r()
      base <- base_r()
      
      mat <- as.matrix(mat)
      storage.mode(mat) <- "integer"
      shiny::validate(
        shiny::need(is.integer(mat), "hf_mat_r() must return an integer matrix"),
        shiny::need(ncol(mat) >= 2L, "shiny::need pop + at least 1 locus."),
        shiny::need(all(mat[,1] > 0, na.rm = TRUE), "Population codes must be positive integers (1..K)"),
        shiny::need(isTRUE(is.finite(base)) && base > 1L, "Invalid base from params."),
        shiny::need(isTRUE(is.finite(conf_level)) && conf_level > 0 && conf_level < 1, "conf_level must be in (0,1)."),
        shiny::need(isTRUE(is.finite(n_boot)) && n_boot >= 10, "n_boot must be >= 10.")
      )
      
      
      # locus names (everything except pop column)
      loci_names <- colnames(mat)[-1L]
      if (is.null(loci_names) || length(loci_names) == 0L) {
        loci_names <- paste0("L", seq_len(ncol(mat) - 1L))
      }
      
      # ---------------------------#
      # 1) Observed (WC84) from C++
      # ---------------------------#
      .bump_progress(15, "Computing observed FST (Weir & Cockerham 1984)...")
      obs_res <- .step("observed_wc84_stats_cpp()", observed_wc84_stats_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      ))
      
      fst_obs_vec <- obs_res$FST
      fst_obs_by_locus <- as.numeric(fst_obs_vec)
      
      loc_obs <- as.character(obs_res$locus_names %||% loci_names)
      names(fst_obs_by_locus) <- loc_obs
      
      # overall FST (WC84 multilocus) = sum(a)/sum(a+b+c)
      fst_obs_overall <- as.numeric(obs_res$FST_overall_ratio_of_sums)
      
      # ---------------------------#
      # 1b) Locus bootstrap (resample loci with replacement)
      # ---------------------------#
      .bump_progress(25, "Bootstrap over loci...")
      locus_comp <- wc84_locus_components_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      )
      # Same algorithm as the former C++ routine, done in R so that the replicate
      # values themselves are available (optional "detail" section of the output).
      locus_boot <- pga_boot_over_loci(
        A = locus_comp$A, B = locus_comp$B, C = locus_comp$C,
        HS = locus_comp$HS, HT = locus_comp$HT,
        n_boot = as.integer(n_boot), conf_level = conf_level, seed = .seed()
      )
      locus_boot_reps <- locus_boot$reps
      locus_boot_tbl  <- data.frame(
        Statistic = locus_boot$table$Statistic,
        Observed  = locus_boot$table$Observed,
        Boot_Mean = colMeans(locus_boot_reps[, locus_boot$table$Statistic, drop = FALSE], na.rm = TRUE),
        SE        = locus_boot$table$SE,
        CI_L      = locus_boot$table$CI_L,
        CI_U      = locus_boot$table$CI_U,
        stringsAsFactors = FALSE
      )

      # ---------------------------#
      # 2) Permutation (single call, fast) from C++
      # ---------------------------#
      perm_res <- NULL
      perm_mat <- NULL
      pvals_by_locus <- rep(NA_real_, length(loci_names))
      pval_overall <- NA_real_
      
      if (!is.null(n_perm) && n_perm > 0) {
        .bump_progress(40, sprintf("Running %s permutations (sub-sample labels shuffled)...",
                                    format(as.integer(n_perm), big.mark = ",")))
        perm_res <-  .step("batch_permute_wc84_fst_auto()", batch_permute_wc84_fst_auto(
          dat            = mat,
          pop_col_1based = 1L,
          missing_code   = as.integer(missing_code),
          base           = as.integer(base),
          B              = as.integer(n_perm),
          n_threads      = .n_threads(),
          seed           = .seed(),
          pval_method    = "greater",              # FST: usually one-sided
          perm_scheme    = "permute_pop_labels",    # breaks pop structure (H0: no structure)
          debug         = FALSE
        ))

        fst_parallel_meta(attr(perm_res, "parallel") %||% NULL)

        perm_mat <- perm_res$FST_perm
        # Ensure colnames align with loci (prefer locus_names from C++ when present)
        if (!is.null(perm_res$locus_names)) colnames(perm_mat) <- as.character(perm_res$locus_names)
        
        pvals_by_locus <- as.numeric(perm_res$p_G)
        if (!is.null(perm_res$locus_names)) names(pvals_by_locus) <- as.character(perm_res$locus_names)
        
        pval_overall <- as.numeric(perm_res$p_FST_overall)
      }
      
      # ---------------------------#
      # 4) Bootstrap (pop blocks) for overall CI (reflected/basic)
      # ---------------------------#
      .bump_progress(70, sprintf("Bootstrap over individuals (%s replicates)...",
                                  format(as.integer(n_boot), big.mark = ",")))
      n_subs_avail <- length(unique(mat[, 1]))
      if (n_subs_avail < 5L) {
        # Fewer than 5 sub-samples (populations): bootstrap over sub-samples
        # is not meaningful - report NA instead of computing.
        nL <- length(loci_names)
        na_mat <- matrix(NA_real_, nrow = as.integer(n_boot), ncol = nL,
                          dimnames = list(NULL, loci_names))
        boot_pop_res <- list(
          locus_names       = loci_names,
          FST_boot          = na_mat, HS_boot = na_mat, HT_boot = na_mat,
          FST_overall_boot  = rep(NA_real_, as.integer(n_boot)),
          HS_overall_boot   = rep(NA_real_, as.integer(n_boot)),
          HT_overall_boot   = rep(NA_real_, as.integer(n_boot)),
          HS_obs            = rep(NA_real_, nL), HT_obs = rep(NA_real_, nL),
          HS_overall_obs    = NA_real_, HT_overall_obs = NA_real_
        )
      } else {
      boot_pop_res <- .step("boot_popblock_wc84_fst_auto()",boot_popblock_wc84_fst_auto(
        mat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base),
        B              = as.integer(n_boot),
        n_threads      = .n_threads(),
        seed           = .seed(),
        debug          = FALSE
      ))
      }

      fst_parallel_meta(attr(boot_pop_res, "parallel") %||% fst_parallel_meta())

      perm_parallel_meta <- if (!is.null(perm_res)) attr(perm_res, "parallel") else NULL
      boot_parallel_meta <- attr(boot_pop_res, "parallel")

      # Pop-block bootstrap matrix (per-locus replicates)
      loc <- as.character(boot_pop_res$locus_names)
      
      boot_pop_mat <- boot_pop_res$FST_boot
      hs_boot_mat  <- boot_pop_res$HS_boot
      ht_boot_mat  <- boot_pop_res$HT_boot
      
      colnames(boot_pop_mat) <- loc
      colnames(hs_boot_mat)  <- loc
      colnames(ht_boot_mat)  <- loc
      
      # Overall bootstrap vector MUST be WC84 multilocus ratio-of-sums.
      boot_overall <- boot_pop_res$FST_overall_boot
      if (is.null(boot_overall) || !is.numeric(boot_overall)) {
        stop("C++ boot_popblock_wc84_fst() did not return numeric FST_overall_boot. Recompile the correct fst_permute_bootstrap.cpp.")
      }
      
      # Summarise using C++ (NORMAL bootstrap CI, FSTAT-like)
      loc <- as.character(boot_pop_res$locus_names)
      
      fst_obs_vec_aligned <- as.numeric(fst_obs_vec)
      names(fst_obs_vec_aligned) <- as.character(obs_res$locus_names %||% loc)
      fst_obs_vec_aligned <- fst_obs_vec_aligned[match(loc, names(fst_obs_vec_aligned))]
      
      sum_pop <- .step("summarize_boot_ci(FST)",summarize_boot_ci(
        boot_mat     = boot_pop_mat,
        obs          = fst_obs_vec_aligned,
        obs_overall  = fst_obs_overall,
        boot_overall = boot_overall,
        confidence   = conf_level
      ))
      
      # =========================#
      # HS / HT CI + means
      # =========================#
      
      # Per-locus bootstrap matrices
      hs_boot_mat <- boot_pop_res$HS_boot
      ht_boot_mat <- boot_pop_res$HT_boot
      
      # Make sure colnames are consistent
      colnames(hs_boot_mat) <- as.character(boot_pop_res$locus_names)
      colnames(ht_boot_mat) <- as.character(boot_pop_res$locus_names)
      
      # Observed vectors (named by locus)
      loc <- as.character(boot_pop_res$locus_names)
      
      hs_obs_vec <- as.numeric(boot_pop_res$HS_obs)
      names(hs_obs_vec) <- loc
      hs_obs_vec <- hs_obs_vec[match(loc, names(hs_obs_vec))]
      
      ht_obs_vec <- as.numeric(boot_pop_res$HT_obs)
      names(ht_obs_vec) <- loc
      ht_obs_vec <- ht_obs_vec[match(loc, names(ht_obs_vec))]
      
      
      # Overall observed + overall bootstrap vectors (from C++)
      hs_obs_overall <- as.numeric(boot_pop_res$HS_overall_obs)
      ht_obs_overall <- as.numeric(boot_pop_res$HT_overall_obs)
      
      hs_overall_boot <- boot_pop_res$HS_overall_boot
      ht_overall_boot <- boot_pop_res$HT_overall_boot
      
      # -- Population-block (subsamples) bootstrap CI --------------------------
      sum_hs <- .step("summarize_boot_ci(HS)", summarize_boot_ci(
        boot_mat     = hs_boot_mat,
        obs          = hs_obs_vec,
        obs_overall  = hs_obs_overall,
        boot_overall = hs_overall_boot,
        confidence   = conf_level
      ))

      sum_ht <- .step("summarize_boot_ci(HT)", summarize_boot_ci(
        boot_mat     = ht_boot_mat,
        obs          = ht_obs_vec,
        obs_overall  = ht_obs_overall,
        boot_overall = ht_overall_boot,
        confidence   = conf_level
      ))

      # -- FST / FIT / FIS / HT : bootstrap over INDIVIDUALS (BS/Indiv) ----------
      # Individuals are resampled with replacement WITHIN each sub-sample (as for
      # FIS); sub-samples are kept. Replaces the former bootstrap over sub-samples
      # for FST, FIT and HT (supervisor's request). FST and HT replicates are
      # biased (resampled genotypes are duplicated), so their percentile CI is
      # shifted by (observed - bootstrap mean), see pga_ci_bias_shift().
      .bump_progress(75, sprintf("Bootstrap over individuals (%s replicates)...",
                                  format(as.integer(n_boot), big.mark = ",")))
      bi <- .step("boot_indiv_wc84_cpp()", boot_indiv_wc84_cpp(
        dat = mat, pop_col_1based = 1L, missing_code = as.integer(missing_code),
        base = as.integer(base), B = as.integer(n_boot), seed = as.numeric(.seed()),
        n_threads = .n_threads()))
      .indiv_sum <- function(boot_mat, obs, obs_overall, boot_overall) {
        colnames(boot_mat) <- loc
        ci <- vapply(seq_along(loc), function(j)
          pga_ci_bias_shift(obs[j], boot_mat[, j], conf_level), numeric(2))
        cio <- pga_ci_bias_shift(obs_overall, boot_overall, conf_level)
        list(mean = colMeans(boot_mat, na.rm = TRUE), ci_lo = ci[1, ], ci_hi = ci[2, ],
             se = apply(boot_mat, 2, stats::sd, na.rm = TRUE),
             overall_mean = mean(boot_overall, na.rm = TRUE),
             overall_ci_lo = cio[["lo"]], overall_ci_hi = cio[["hi"]],
             overall_se = stats::sd(boot_overall, na.rm = TRUE))
      }
      ht_obs_vec2 <- as.numeric(obs_res$HT)[match(loc, as.character(obs_res$locus_names %||% loc))]
      ht_obs_ov2  <- if (is.finite(ht_obs_overall)) ht_obs_overall else mean(ht_obs_vec2, na.rm = TRUE)
      sum_pop <- .indiv_sum(bi$FST_boot, fst_obs_vec_aligned, fst_obs_overall, bi$FST_overall_boot)
      sum_ht  <- .indiv_sum(bi$HT_boot,  ht_obs_vec2,          ht_obs_ov2,     bi$HT_overall_boot)
      ht_obs_vec <- ht_obs_vec2; ht_obs_overall <- ht_obs_ov2
      ht_boot_mat <- bi$HT_boot; colnames(ht_boot_mat) <- loc
      ht_overall_boot <- bi$HT_overall_boot
      boot_pop_mat <- bi$FST_boot; colnames(boot_pop_mat) <- loc
      boot_overall <- bi$FST_overall_boot
      fis_obs_ov <- as.numeric(obs_res$FIS_overall_ratio_of_sums)
      fit_obs_ov <- as.numeric(obs_res$FIT_overall_ratio_of_sums)
      ci_fit_ov <- pga_ci_bias_shift(fit_obs_ov, bi$FIT_overall_boot, conf_level)
      ci_fis_ov <- pga_ci_bias_shift(fis_obs_ov, bi$FIS_overall_boot, conf_level)
      subs_fitfis <- list(FIT = ci_fit_ov, FIS = ci_fis_ov,
                          reps = cbind(FIT = bi$FIT_overall_boot, FIS = bi$FIS_overall_boot))

      # -- Individual bootstrap CI for HS (resample individuals within pops) --
      # Genetic Diversities only (skipped for Subdivision, see diversity_extras).
      if (isTRUE(diversity_extras)) {
        indiv_boot_res <- .step("boot_indiv_hs_cpp()", boot_indiv_hs_cpp(
          dat            = mat,
          pop_col_1based = 1L,
          missing_code   = as.integer(missing_code),
          base           = as.integer(base),
          B              = as.integer(n_boot),
          seed           = .seed(),
          n_threads      = .n_threads()
        ))
        hs_indiv_boot_mat <- indiv_boot_res$HS_boot
        colnames(hs_indiv_boot_mat) <- as.character(indiv_boot_res$locus_names)
        hs_indiv_overall_boot <- indiv_boot_res$HS_overall_boot

        sum_hs_indiv <- .step("summarize_boot_ci(HS indiv)", summarize_boot_ci(
          boot_mat     = hs_indiv_boot_mat,
          obs          = hs_obs_vec,
          obs_overall  = hs_obs_overall,
          boot_overall = hs_indiv_overall_boot,
          confidence   = conf_level
        ))
      }

      loc_fst <- loc
      loc_hs  <- loc
      loc_ht  <- loc

      if (!is.null(names(fst_obs_by_locus))) {
        fst_obs_by_locus <- fst_obs_by_locus[match(loc_fst, names(fst_obs_by_locus))]
        names(fst_obs_by_locus) <- loc_fst
      }
      if (!is.null(names(pvals_by_locus))) {
        pvals_by_locus <- pvals_by_locus[match(loc_fst, names(pvals_by_locus))]
        names(pvals_by_locus) <- loc_fst
      }

      # -- Helper: build per-locus + Overall HS table for one bootstrap mode --
      # bias_shift = TRUE (bootstrap over INDIVIDUALS): the percentile CI is shifted by
      # (observed - bootstrap mean), see pga_ci_bias_shift() for the reason.
      .hs_boot_tbl <- function(sum_obj, indiv_boot_mat, indiv_boot_overall,
                               obs_vec, obs_overall, loc_names, bias_shift = FALSE) {
        boot_se_loci   <- apply(indiv_boot_mat, 2, sd, na.rm = TRUE)
        boot_se_overall <- sd(as.numeric(indiv_boot_overall), na.rm = TRUE)
        per_locus <- data.frame(
          ID          = loc_names,
          Observed_HS = as.numeric(obs_vec[loc_names]),
          Boot_Mean   = as.numeric(sum_obj$mean),
          Boot_SE     = as.numeric(boot_se_loci),
          CI_L        = as.numeric(sum_obj$ci_lo),
          CI_U        = as.numeric(sum_obj$ci_hi),
          stringsAsFactors = FALSE
        )
        if (isTRUE(bias_shift)) {
          sh <- per_locus$Observed_HS - per_locus$Boot_Mean
          per_locus$CI_L <- per_locus$CI_L + sh
          per_locus$CI_U <- per_locus$CI_U + sh
        }
        overall <- data.frame(
          ID          = "Overall",
          Observed_HS = as.numeric(obs_overall),
          Boot_Mean   = as.numeric(sum_obj$overall_mean),
          Boot_SE     = as.numeric(boot_se_overall),
          CI_L        = as.numeric(sum_obj$overall_ci_lo),
          CI_U        = as.numeric(sum_obj$overall_ci_hi),
          stringsAsFactors = FALSE
        )
        if (isTRUE(bias_shift)) {
          sh <- overall$Observed_HS - overall$Boot_Mean
          overall$CI_L <- overall$CI_L + sh
          overall$CI_U <- overall$CI_U + sh
        }
        rbind(per_locus, overall)
      }

      # Table 1: HS by individuals (Genetic Diversities only)
      hs_indiv_tbl <- if (isTRUE(diversity_extras))
        .hs_boot_tbl(sum_hs_indiv,
                     hs_indiv_boot_mat, hs_indiv_overall_boot,
                     hs_obs_vec, hs_obs_overall, loc_hs, bias_shift = TRUE)
      else NULL

      # Table 2: HS by populations (block bootstrap)
      hs_pop_tbl <- .hs_boot_tbl(sum_hs,
                                  hs_boot_mat, hs_overall_boot,
                                  hs_obs_vec, hs_obs_overall, loc_hs)

      # Table 3: HS by loci - Overall only (from locus_boot_tbl)
      hs_locus_tbl <- tryCatch({
        row <- locus_boot_tbl[locus_boot_tbl$Statistic == "HS", ,drop = FALSE]
        data.frame(
          Statistic   = "HS",
          Observed    = as.numeric(row$Observed),
          Boot_Mean   = as.numeric(row$Boot_Mean),
          Boot_SE     = as.numeric(row$SE),
          CI_L        = as.numeric(row$CI_L),
          CI_U        = as.numeric(row$CI_U),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(Statistic="HS", Observed=hs_obs_overall,
                   Boot_Mean=NA_real_, Boot_SE=NA_real_,
                   CI_L=NA_real_, CI_U=NA_real_,
                   stringsAsFactors=FALSE)
      })

      # HT tables (populations + locus bootstrap, unchanged logic)
      ht_pop_se_loci   <- apply(ht_boot_mat, 2, sd, na.rm = TRUE)
      ht_pop_se_overall <- sd(as.numeric(ht_overall_boot), na.rm = TRUE)
      ht_tbl <- data.frame(
        ID          = loc_ht,
        Observed_HT = as.numeric(ht_obs_vec[loc_ht]),
        Boot_Mean   = as.numeric(sum_ht$mean),
        Boot_SE     = as.numeric(ht_pop_se_loci),
        CI_L        = as.numeric(sum_ht$ci_lo),
        CI_U        = as.numeric(sum_ht$ci_hi),
        stringsAsFactors = FALSE
      )
      ht_overall_row <- data.frame(
        ID          = "Overall",
        Observed_HT = as.numeric(ht_obs_overall),
        Boot_Mean   = as.numeric(sum_ht$overall_mean),
        Boot_SE     = as.numeric(ht_pop_se_overall),
        CI_L        = as.numeric(sum_ht$overall_ci_lo),
        CI_U        = as.numeric(sum_ht$overall_ci_hi),
        stringsAsFactors = FALSE
      )
      ht_final <- rbind(ht_tbl, ht_overall_row)

      # Table 4: HS per population (Genetic Diversities only; skipped for Subdivision).
      # Observed HS = mean over loci weighted by the number of individuals typed at each
      # locus (Genepop convention), bootstrap over individuals AND over loci:
      # see pga_hs_per_population() in utils_stats.R.
      hs_per_pop_res <- NULL
      hs_per_pop_tbl <- NULL
      if (isTRUE(diversity_extras)) {
        hs_per_pop_res <- pga_hs_per_population(mat, base = base, n_boot = n_boot,
                                                conf_level = conf_level, seed = .seed(),
                                                missing_code = missing_code)
        hs_per_pop_tbl <- hs_per_pop_res$table
      }

      # Backward-compat alias used by value box and download handler
      hs_final <- hs_pop_tbl
      
      boot_mean <- as.numeric(sum_pop$mean)
      ci_lower  <- as.numeric(sum_pop$ci_lo)
      ci_upper  <- as.numeric(sum_pop$ci_hi)
      
      names(boot_mean) <- loc_fst
      names(ci_lower)  <- loc_fst
      names(ci_upper)  <- loc_fst
      
      # Overall outputs
      overall_boot_mean <- as.numeric(sum_pop$overall_mean)
      overall_ci_l      <- as.numeric(sum_pop$overall_ci_lo)
      overall_ci_u      <- as.numeric(sum_pop$overall_ci_hi)
      
      # Optional: median from the bootstrap matrix (if you still want it)
      boot_median <- apply(boot_pop_mat, 2, median, na.rm = TRUE)
      
      
      # ---------------------------#
      # 6) Final results table
      # ---------------------------#
      loc <- loc_fst
      if (is.null(loc) || length(loc) == 0L) loc <- loci_names
      
      obs_vec <- fst_obs_by_locus
      if (!is.null(names(obs_vec))) obs_vec <- obs_vec[match(loc, names(obs_vec))]
      
      p_vec <- pvals_by_locus
      if (!is.null(names(p_vec))) p_vec <- p_vec[match(loc, names(p_vec))]
      
      res_tbl <- data.frame(
        ID          = loc,
        Observed_FST = as.numeric(obs_vec),
        Boot_Mean    = as.numeric(boot_mean[loc]),
        Boot_Median  = as.numeric(boot_median[loc]),
        SE           = as.numeric(sum_pop$se[loc]),
        P_value      = as.numeric(p_vec),
        CI_L         = as.numeric(ci_lower[loc]),
        CI_U         = as.numeric(ci_upper[loc]),
        stringsAsFactors = FALSE
      )
      
      # overall_boot_mean already computed from v_overall above (more coherent)
      # if it was NA (too few finite), fallback:
      # Option 1: define an overall "bootstrap median" coherently from boot_overall
      overall_boot_median <- if (exists("boot_overall") && length(boot_overall) > 0) {
        median(boot_overall, na.rm = TRUE)
      } else {
        NA_real_
      }
      
      overall_row <- data.frame(
        ID          = "Overall",
        Observed_FST = fst_obs_overall,
        Boot_Mean    = overall_boot_mean,
        Boot_Median  = overall_boot_median,
        SE           = sum_pop$overall_se,
        P_value      = pval_overall,
        CI_L         = overall_ci_l,
        CI_U         = overall_ci_u,
        stringsAsFactors = FALSE
      )
      
      # force same column order as res_tbl
      overall_row <- overall_row[, names(res_tbl), drop = FALSE]
      
      final_results <- rbind(res_tbl, overall_row)

      # W&C theta permutation test (individuals reassigned among sub-samples):
      # P_theta_gt = P(permuted theta >= observed)  [H1: theta > 0, structure]
      # P_theta_lt = P(permuted theta <= observed)  [H1: theta < 0]
      final_results$P_theta_gt <- NA_real_; final_results$P_theta_lt <- NA_real_
      if (!is.null(perm_res) && !is.null(perm_res$FST_perm)) {
        .pp <- function(perm, obs) {
          perm <- perm[is.finite(perm)]
          if (!is.finite(obs) || !length(perm)) return(c(NA_real_, NA_real_))
          c((pga_n_ge(perm, obs) + 1) / (length(perm) + 1), (pga_n_le(perm, obs) + 1) / (length(perm) + 1))
        }
        pm <- perm_res$FST_perm
        for (i in seq_along(loc)) {
          j <- match(loc[i], colnames(pm)); if (is.na(j)) j <- i
          v <- .pp(pm[, j], final_results$Observed_FST[i])
          final_results$P_theta_gt[i] <- v[1]; final_results$P_theta_lt[i] <- v[2]
        }
        v <- .pp(as.numeric(perm_res$FST_overall_perm), fst_obs_overall)
        n_f <- nrow(final_results)
        final_results$P_theta_gt[n_f] <- v[1]; final_results$P_theta_lt[n_f] <- v[2]
      }

      .bump_progress(95, "Finalising results...")

      list(
        final_table      = final_results,
        locus_boot_table = locus_boot_tbl,
        # Overall (multilocus) CI by bootstrap over sub-samples, to sit next to the
        # bootstrap-over-loci CI in the multilocus table.
        subs_overall = {
          hs_ci <- pga_ci_bias_shift(hs_obs_overall, bi$HS_overall_boot, conf_level)
          data.frame(
            Statistic = c("FIS", "FST", "FIT", "HS", "HT"),
            SE_indiv  = c(stats::sd(bi$FIS_overall_boot), sum_pop$overall_se, stats::sd(bi$FIT_overall_boot),
                          stats::sd(bi$HS_overall_boot), sum_ht$overall_se),
            CI_L_subs = c(subs_fitfis$FIS[["lo"]], sum_pop$overall_ci_lo, subs_fitfis$FIT[["lo"]],
                          hs_ci[["lo"]], sum_ht$overall_ci_lo),
            CI_U_subs = c(subs_fitfis$FIS[["hi"]], sum_pop$overall_ci_hi, subs_fitfis$FIT[["hi"]],
                          hs_ci[["hi"]], sum_ht$overall_ci_hi),
            stringsAsFactors = FALSE)
        },
        # Replicate values of every bootstrap (optional "detail" sections of the files).
        boot_detail = list(
          fst_subs     = cbind(boot_pop_mat, Overall = boot_overall),
          hs_subs      = cbind(hs_boot_mat,  Overall = hs_overall_boot),
          ht_subs      = cbind(ht_boot_mat,  Overall = ht_overall_boot),
          hs_indiv     = if (isTRUE(diversity_extras)) cbind(hs_indiv_boot_mat, Overall = hs_indiv_overall_boot) else NULL,
          hs_pop_indiv = if (!is.null(hs_per_pop_res)) hs_per_pop_res$reps_indiv else NULL,
          hs_pop_loci  = if (!is.null(hs_per_pop_res)) hs_per_pop_res$reps_loci  else NULL,
          loci         = locus_boot_reps,
          fit_fis_subs = subs_fitfis$reps
        ),
        hs_table     = hs_final,
        hs_indiv_tbl = hs_indiv_tbl,
        hs_pop_tbl   = hs_pop_tbl,
        hs_locus_tbl = hs_locus_tbl,
        hs_per_pop_tbl = hs_per_pop_tbl,
        ht_table    = ht_final,
        observed_fst = fst_obs_by_locus,
        permutation_results = perm_mat,
        boot_parallel_meta = boot_parallel_meta,
        perm_parallel_meta = perm_parallel_meta,
        bootstrap_results = list(
          population = boot_pop_res,
          population_matrix = boot_pop_mat,
          overall_boot = boot_overall,
          hs = sum_hs,
          ht = sum_ht
        ),
        metadata = list(
          parallel_perm = perm_parallel_meta,
          parallel_boot = boot_parallel_meta,
          requested_threads = as.integer(.n_threads()),
          n_loci = length(loci_names),
          n_permutations = n_perm,
          n_bootstrap = n_boot,
          conf_level = conf_level,
          base = base,
          loci_names = loci_names,
          pop_names = as.character(attr(mat, "pop_levels")),
          dataset_name = if (!is.null(rv$dataset_filename)) rv$dataset_filename else NA_character_,
          pval_method = if (!is.null(perm_res)) perm_res$pval_method else NA_character_
        )
        
      )
    }
    
    
    ## Run FST bootstrap and permutation (button) ----
    .run_subdivision_fst_computation <- function() {
      
      db_ready()
      
      if (input$n_boot_fst < 10) {
        showNotification("Number of bootstrap replicates should be at least 10", type = "warning")
        return(NULL)
      }
      
      results <- tryCatch({
        start_time <- Sys.time()
        
        shinyWidgets::updateProgressBar(session, "fst_progress", value = 0,
                                         title = "Starting FST analysis...")
        
        results <- run_bootstrap_fst_analysis(
          n_perm         = input$n_perm_fst,
          n_boot         = input$n_boot_fst,
          conf_level     = input$conf_level_fst,
          missing_code   = 0L,
          progress_id    = "fst_progress",
          diversity_extras = FALSE      # HS-by-individuals / HS-per-population: Diversities only
        )
        
        shinyWidgets::updateProgressBar(session, "fst_progress", value = 100,
                                         title = "Done")
        
        duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
        fst_boot_timing(duration)
        fst_boot_results(results)
        
        showNotification(
          "Computations completed.",
          type = "message"
        )
        results
        
      }, error = function(e) {
        fst_boot_results(NULL)
        fst_boot_timing(NULL)
        showNotification(paste("Error in FST analysis:", e$message), type = "error")
        NULL
      })

      results
    }

    ## Run FST bootstrap and permutation (button from Genetic diversities tab) ----
    .run_diversities_computation <- function() {

      db_ready()

      if (input$n_boot_fst_div < 10) {
        showNotification("Number of bootstrap replicates should be at least 10", type = "warning")
        return(NULL)
      }

      results <- tryCatch({
        start_time <- Sys.time()

        shinyWidgets::updateProgressBar(session, "fst_progress_div", value = 0,
                                         title = "Starting FST analysis...")

        results <- run_bootstrap_fst_analysis(
          n_perm         = 0L,   # no permutation test in Diversities
          n_boot         = input$n_boot_fst_div,
          conf_level     = 0.95,   # fixed: the Diversities module has no confidence-level widget
          missing_code   = 0L,
          progress_id    = "fst_progress_div"
        )

        shinyWidgets::updateProgressBar(session, "fst_progress_div", value = 100,
                                         title = "Done")

        duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
        div_boot_timing(duration)
        div_boot_results(results)

        showNotification(
          "Computations completed.",
          type = "message"
        )
        results

      }, error = function(e) {
        div_boot_results(NULL)
        div_boot_timing(NULL)
        showNotification(paste("Error in FST analysis:", e$message), type = "error")
        NULL
      })

      results
    }

    ## ===== FST value boxes =====
    



    # -- Methods block written at the TOP of the single Subdivision output file.
    .write_fst_params <- function(con, res, gres) {
      md <- if (is.list(res)) res$metadata else NULL
      loci <- md$loci_names %||% character(0)
      pops <- md$pop_names  %||% character(0)
      nb   <- md$n_bootstrap %||% input$n_boot_fst
      npm  <- md$n_permutations %||% input$n_perm_fst
      nc   <- if (is.list(gres)) gres$metadata$n_complete else NULL
      hdr <- c(
        "Population Subdivision - FST (Weir & Cockerham 1984) and G-based test (Goudet et al. 1996)",
        sprintf("Dataset: %s", if (!is.null(md$dataset_name) && !is.na(md$dataset_name)) md$dataset_name else "default_dataset"),
        sprintf("Loci (n = %d): %s", length(loci), paste(loci, collapse = ", ")),
        sprintf("Sub-samples (n = %d): %s", length(pops), paste(pops, collapse = ", ")),
        "",
        sprintf("P-values of Section 1: one-sided permutation test based on the log-likelihood G statistic (as in FSTAT, Goudet et al. 1996), %s permutations; the individuals typed at a locus are reassigned at random among sub-samples (sub-sample labels shuffled).", npm),
        "  p = (b + 1) / (m + 1), b = number of permuted G >= observed G (ties included), m = number of permutations; the Overall p-value uses G summed over loci.",
        "",
        sprintf("Confidence intervals: %s level, percentile bootstrap, %s replicates for each type of bootstrap:", md$conf_level %||% input$conf_level_fst, nb),
        "  BS/Indiv, _indiv = bootstrap over INDIVIDUALS (individuals resampled with replacement within each sub-sample; sub-samples kept). Replicates of FST and HT are slightly biased (duplicated genotypes), so the percentile interval is shifted by (observed - bootstrap mean); its width is the bootstrap one.",
        "  BS/Loci, _loci = bootstrap over LOCI (loci resampled with replacement).",
        "  SE = bootstrap standard error (standard deviation of the bootstrap replicates).",
        "  P_theta_gt / P_theta_lt (Section 1): permutation test on the W&C estimator theta (FST): proportion of permuted theta >= observed (H1: theta > 0) and <= observed (H1: theta < 0), (b + 1)/(m + 1).",
        "",
        sprintf("G-based test (Section 2): %s permutations of COMPLETE MULTILOCUS GENOTYPES (whole individuals) among sub-samples; valid when Hardy-Weinberg is NOT assumed within samples.",
                if (is.list(gres)) gres$metadata$n_perm else npm),
        "  p_>= = proportion of permuted G >= observed G; p_> = proportion of permuted G > observed G (FSTAT gives both); (b + 1) / (m + 1).",
        if (!is.null(nc)) sprintf("  Number of complete multilocus genotypes used, per sample: %s; total = %d",
                                   paste(sprintf("%s = %d", names(nc), nc), collapse = ", "), sum(nc)),
        "  Overall row = G summed over loci, tested the same way.",
        ""
      )
      writeLines(hdr, con = con, useBytes = TRUE)
    }

    # -- One button, one click, one action: clicking "Run" IS the download
    #    request itself - the FST bootstrap/permutation AND the G-based test
    #    run inside this same content() function, and everything (methods
    #    block + FST per locus + multilocus estimators + G-based test, then the
    #    optional bootstrap replicate values) is written to ONE plain .txt
    #    file (no zip, no figure).
    #    Populates fst_boot_results()/fst_boot_timing() - Subdivision's own
    #    value boxes only; Genetic Diversities has its own separate
    #    div_boot_results()/div_boot_timing().
    output$ui_fst_out_status <- renderUI({
      tags$p(style = "color:#555;font-size:14px;margin-top:6px;",
        "The results will be saved in ", tags$code(paste0("subdivision_FST_", Sys.Date(), ".txt")), ".")
    })

    output$run_FST_Analysis <- downloadHandler(
      filename = function() paste0("subdivision_FST_", Sys.Date(), ".txt"),
      content  = function(file) {
        res <- .run_subdivision_fst_computation()
        req(res)
        # G-based test: same number of permutations as the FST test.
        gres <- .run_g_test_computation(n_perm = input$n_perm_fst)
        nb   <- res$metadata$n_bootstrap
        npm  <- res$metadata$n_permutations

        con <- file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        .write_fst_params(con, res, gres)

        writeLines(sprintf("Section 1: FST per locus - confidence intervals with %s bootstrap over individuals (locus rows and Overall BS/Indiv) and over loci (Overall BS/Loci); one-sided permutation p-values based on the G statistic (%s permutations)", nb, npm), con = con)
        pga_write_table(pga_fst_section(res$final_table, res$locus_boot_table), file = con)
        writeLines("", con = con)

        if (is.list(gres)) {
          gt <- gres$final_table[, c("ID", "G_obs", "p_ge", "p_gt")]
          names(gt) <- c("ID", "G_obs", "p_>=", "p_>")
          writeLines(sprintf("Section 2: G-based permutation test of subdivision (multilocus genotypes permuted among sub-samples), %s permutations", gres$metadata$n_perm), con = con)
          pga_write_table(gt, file = con)
        } else {
          writeLines("Section 2: G-based test could not be computed (see the notification shown in the app).", con = con)
        }

        # Pairwise FST between sub-samples: no GPS needed, always produced.
        pw <- tryCatch(pga_pairwise_fst(hf_mat_r(), base = base_r(), n_boot = nb,
                                        n_perm = npm, conf_level = input$conf_level_fst,
                                        seed = .seed(), n_threads = .n_threads()),
                       error = function(e) { showNotification(paste("Pairwise FST failed:", conditionMessage(e)), type = "warning"); NULL })
        writeLines("", con = con)
        if (!is.null(pw)) {
          writeLines(sprintf("Section 3: pairwise FST between sub-samples (multilocus) - CI with %s bootstrap over individuals (bias-shifted) and over loci; one-sided permutation p-values (%s permutations of individuals between the two sub-samples)", nb, npm), con = con)
          pga_write_table(pw$overall, file = con)
          writeLines("", con = con)
          writeLines(sprintf("Section 4: pairwise FST between sub-samples, per locus - CI with %s bootstrap over individuals (bias-shifted); one-sided permutation p-values (G statistic)", nb), con = con)
          pga_write_table(pw$per_locus, file = con)
        }

        if (isTRUE(input$fst_detail)) {
          writeLines("", con = con)
          writeLines(sprintf("DETAIL: the %s bootstrap replicate values behind the confidence intervals above", nb), con = con)
          writeLines("", con = con)
          pga_write_replicates(con, "Detail A: FST, bootstrap over individuals (raw values before the bias shift; one column per locus, then Overall)", res$boot_detail$fst_subs)
          pga_write_replicates(con, "Detail B: multilocus FST, FIT, FIS, HS, HT, bootstrap over loci", res$boot_detail$loci)
        }
      }
    )






    
    
    ## ===== Diversities tab value boxes (global_fst, global_hs, global_ht, time) =====






    



    

    
    
    ## ===== FST, HT, HS  plots =====

    # Helper: build a combined per-locus + Overall plot for HS or HT.
    # Mirrors FST plot style: size=3, width=0.2, no size/linewidth aesthetics.

    # -- Methods block written at the TOP of the single Diversities file.
    .write_div_params <- function(con, res) {
      md   <- if (is.list(res)) res$metadata else NULL
      loci <- md$loci_names %||% character(0)
      pops <- md$pop_names  %||% character(0)
      nb   <- md$n_bootstrap %||% input$n_boot_fst_div
      hdr <- c(
        "Diversities confidence intervals - HS and HT (Nei & Chesser 1983)",
        sprintf("Dataset: %s", if (!is.null(md$dataset_name) && !is.na(md$dataset_name)) md$dataset_name else "default_dataset"),
        sprintf("Loci (n = %d): %s", length(loci), paste(loci, collapse = ", ")),
        sprintf("Sub-samples (n = %d): %s", length(pops), paste(pops, collapse = ", ")),
        "",
        "Estimators: HS = unbiased gene diversity within populations, n/(n-1) * (1 - sum(p^2) - Ho/(2n)) per population and locus, averaged over populations;",
        "            HT = total gene diversity from the mean allele frequencies.",
        "            HS per sub-sample: mean over loci weighted by the number of individuals typed at each locus (as in Genepop).",
        "",
        sprintf("Confidence intervals: %s level, percentile bootstrap, %s replicates for each type of bootstrap:", md$conf_level %||% 0.95, nb),
        "  _indiv = bootstrap over INDIVIDUALS (individuals resampled within each sub-sample).",
        "           Resampling individuals makes every replicate slightly less diverse than the sample (the bootstrap mean is about HS/(2n) below the observed HS),",
        "           so these intervals are shifted by (observed HS - bootstrap mean); their width is the bootstrap one.",
        "  _subs  = bootstrap over SUB-SAMPLES (populations resampled as blocks); NA if fewer than 5 sub-samples. Used for HS per locus only.",
        "  _loci  = bootstrap over LOCI (loci resampled with replacement).",
        "  SE = bootstrap standard error (standard deviation of the bootstrap replicates).",
        "  Multilocus HS and HT (Section 5, _indiv): recomputed on every replicate.",
        ""
      )
      writeLines(hdr, con = con, useBytes = TRUE)
    }

    # Writes one titled section: a title line, then the table (no `append`
    # argument: writing to an already-open connection just continues).
    .write_div_section <- function(con, title, df) {
      writeLines(title, con = con, useBytes = TRUE)
      pga_write_table(df, file = con)
      writeLines("", con = con)
    }

    # -- One button, one click, one action: clicking "Run" IS the download
    #    request itself - the bootstrap analysis runs inside this same
    #    content() function and ALL results are written to ONE plain .txt
    #    file (methods block + 5 sections + optional replicate values; no zip,
    #    no figures).
    output$ui_div_out_status <- renderUI({
      tags$p(style = "color:#555;font-size:14px;margin-top:6px;",
        "The results will be saved in ", tags$code(paste0("Diversities_confidence_intervals_", Sys.Date(), ".txt")), ".")
    })

    output$run_FST_Analysis_div <- downloadHandler(
      filename = function() paste0("Diversities_confidence_intervals_", Sys.Date(), ".txt"),
      content  = function(file) {
        res <- .run_diversities_computation()
        req(res)
        nb <- res$metadata$n_bootstrap

        con <- file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        .write_div_params(con, res)
        .write_div_section(con, sprintf("Section 1: HS per locus - confidence intervals with %s bootstrap over individuals (resampled within each sub-sample, bias-shifted)", nb),
                           pga_diversity_section(res$hs_indiv_tbl, "Observed_HS"))
        .write_div_section(con, sprintf("Section 2: HS per locus - confidence intervals with %s bootstrap over sub-samples (sub-samples resampled as blocks)", nb),
                           pga_diversity_section(res$hs_pop_tbl, "Observed_HS"))
        .write_div_section(con, sprintf("Section 3: HS per sub-sample - confidence intervals with %s bootstrap over individuals (bias-shifted) and over loci; Overall row: HS over all loci and sub-samples", nb),
                           pga_hs_pop_section(res$hs_per_pop_tbl, res$hs_indiv_tbl, res$locus_boot_table))
        .write_div_section(con, sprintf("Section 4: HT per locus - confidence intervals with %s bootstrap over individuals (bias-shifted)", nb),
                           pga_diversity_section(res$ht_table, "Observed_HT"))
        .write_div_section(con, sprintf("Section 5: HS and HT (multilocus) - confidence intervals with %s bootstrap over loci and over individuals", nb),
                           pga_multilocus_section(res$locus_boot_table, res$subs_overall, stats = c("HS", "HT")))

        if (isTRUE(input$div_detail)) {
          writeLines(sprintf("DETAIL: the %s bootstrap replicate values behind the confidence intervals above", nb), con = con)
          writeLines("", con = con)
          bd <- res$boot_detail
          pga_write_replicates(con, "Detail A: HS per locus, bootstrap over individuals (one column per locus, then Overall; raw values, before the bias shift)", bd$hs_indiv)
          pga_write_replicates(con, "Detail B: HS per locus, bootstrap over sub-samples", bd$hs_subs)
          pga_write_replicates(con, "Detail C: HS per sub-sample, bootstrap over individuals (raw values, before the bias shift)", bd$hs_pop_indiv)
          pga_write_replicates(con, "Detail D: HS per sub-sample, bootstrap over loci", bd$hs_pop_loci)
          pga_write_replicates(con, "Detail E: HT per locus, bootstrap over individuals (raw values before the bias shift)", bd$ht_subs)
          pga_write_replicates(con, "Detail F: multilocus FST, FIT, FIS, HS, HT, bootstrap over loci", bd$loci)
        }
      }
    )




    ## Reactive containers ----

    ## G a partir d'un tableau de comptage n_pop x n_allele ("flat", vecteur) ----
    .g_stat_from_flat <- function(cnt_flat, n_pop, n_allele) {
      cnt <- matrix(cnt_flat, nrow = n_pop, ncol = n_allele)
      rs  <- rowSums(cnt); cs <- colSums(cnt); n <- sum(cnt)
      if (n == 0L || sum(rs > 0L) < 2L || sum(cs > 0L) < 2L) return(NA_real_)
      E  <- outer(rs, cs) / n
      ok <- cnt > 0L & E > 0
      if (!any(ok)) return(NA_real_)
      2 * sum(cnt[ok] * log(cnt[ok] / E[ok]))
    }

    ## Tableau allele x population vectorise via tabulate() ----
    ## pop_idx0    : code population 0-based, longueur n_valid
    ## allele_idx0 : code allele 0-based pour c(a1,a2), longueur 2*n_valid (PRE-CALCULE, fixe)
    .g_stat_vec <- function(pop_idx0, allele_idx0, n_pop, n_allele) {
      pop2 <- c(pop_idx0, pop_idx0)                       # meme pop pour a1 et a2 d'un individu
      bin  <- pop2 + n_pop * allele_idx0 + 1L              # index 1-based dans la matrice n_pop x n_allele
      cnt_flat <- tabulate(bin, nbins = n_pop * n_allele)  # comptage en C - rapide
      .g_stat_from_flat(cnt_flat, n_pop, n_allele)
    }

    ## Observer: Run button ----
    # Runs the FSTAT-style G-based subdivision test and RETURNS the result list
    # (NULL on failure). Called from the Subdivision Run button, with the same
    # number of permutations as the FST test (no separate G-test parameters).
    .run_g_test_computation <- function(n_perm) {
      db_ready()
      n_perm <- as.integer(n_perm)
      if (!is.finite(n_perm) || n_perm < 1000L) {
        showNotification("Minimum 1 000 permutations required for the G-based test.", type = "warning")
        return(NULL)
      }

      ok <- tryCatch({
        start_time <- Sys.time()

        # -- Sources DB-first (meme pattern que FST) --------------------------
        mat  <- hf_mat_r()   # colonnes deja reordonnees via loci_order_r() dans hf_mat_r
        base <- base_r()

        # NOTE: hf_mat_r() already returns a genuine matrix carrying a
        # "pop_levels" attribute (the real population names), so no
        # as.matrix() is needed here; storage.mode<- below keeps attributes
        # and is enough to guarantee integer type. (The generic "Pop1".."PopN"
        # labels once seen in the G-test output came from the column
        # re-ordering inside hf_mat_r(), which dropped the attribute - fixed
        # there, see the comment in hf_mat_r().)
        storage.mode(mat) <- "integer"
        shiny::validate(
          shiny::need(is.integer(mat),                      "hf_mat_r() must return an integer matrix"),
          shiny::need(ncol(mat) >= 2L,                      "Need pop + at least 1 locus"),
          shiny::need(all(mat[, 1L] > 0, na.rm = TRUE),     "Population codes must be positive integers"),
          shiny::need(isTRUE(is.finite(base)) && base > 1L, "Invalid base from params")
        )

        # loci_names dans l'ordre physique DuckDB (hf_mat_r() reordonne par loci_order_r())
        loci_names <- colnames(mat)[-1L]
        if (is.null(loci_names) || length(loci_names) == 0L)
          loci_names <- paste0("L", seq_len(ncol(mat) - 1L))

        n_loci <- length(loci_names)

        pop_codes <- as.integer(mat[, 1L])

        # -- Restriction FSTAT : genotypes multi-locus COMPLETS uniquement --------
        # FSTAT (Goudet et al. 1996, ?7.1, note 1) : "Only complete multilocus
        # genotypes are randomised, to make sure that the combined test over loci
        # is valid." L'en-tete des fichiers FSTAT_G ("Number of complete
        # multilocus genotypes in the different samples") confirme qu'un seul
        # sous-ensemble FIXE d'individus (non manquants a TOUS les loci
        # simultanement) est utilise pour l'ensemble du test, et pas,
        # comme precedemment ici, un sous-ensemble different par locus selon
        # ses propres donnees manquantes (cela change a la fois le G observe
        # et la distribution nulle de permutation).
        geno_mat      <- mat[, -1L, drop = FALSE]
        complete_mask <- rowSums(is.na(geno_mat) | geno_mat <= 0L) == 0L
        shiny::validate(shiny::need(
          sum(complete_mask) >= 2L,
          "Not enough individuals with a complete multilocus genotype across all loci."
        ))

        pop_codes_complete <- pop_codes[complete_mask]
        pops      <- sort(unique(pop_codes_complete[is.finite(pop_codes_complete) & pop_codes_complete > 0L]))
        n_pops    <- length(pops)
        pop_idx0_full <- match(pop_codes_complete, pops) - 1L   # 0-based, sous-ensemble complet
        mat_complete  <- mat[complete_mask, , drop = FALSE]

        # Noms de populations : source UNIQUE (attribut pose par hf_mat_r()) - plus de
        # requete DB separee dont l'ordre alphabetique n'est pas garanti correspondre
        # positionnellement aux codes `pops` (source du mauvais etiquetage precedent).
        pop_levels_attr <- attr(mat, "pop_levels")
        pop_names <- if (!is.null(pop_levels_attr)) as.character(pop_levels_attr)[pops]
                     else paste0("Pop", pops)

        # -- Decoder les genotypes et pre-calculer les index alleles (fixes, hors permutation) --
        # gt = a1*base + a2 - identique au decodage dans ld_pvalues_cpp
        # Tous les individus restants (mat_complete) sont, par construction,
        # non manquants a CHAQUE locus - loci_idx est donc identique (1:n) pour
        # tous les loci, mais on garde le meme schema generique par securite.
        loci_idx      <- vector("list", n_loci)  # lignes de mat_complete valides pour ce locus
        loci_allele0  <- vector("list", n_loci)  # code allele 0-based, longueur 2*n_valid (fixe)
        loci_n_allele <- integer(n_loci)
        loci_n_valid  <- integer(n_loci)

        for (j in seq_len(n_loci)) {
          g   <- as.integer(mat_complete[, j + 1L])
          ok  <- is.finite(g) & g > 0L
          a1  <- g[ok] %/% base
          a2  <- g[ok] %% base
          ok2 <- a1 > 0L & a2 > 0L
          a1  <- a1[ok2]; a2 <- a2[ok2]

          idx      <- which(ok)[ok2]
          alleles  <- sort(unique(c(a1, a2)))
          n_allele <- length(alleles)

          loci_idx[[j]]     <- idx
          loci_allele0[[j]] <- match(c(a1, a2), alleles) - 1L
          loci_n_allele[j]  <- n_allele
          loci_n_valid[j]   <- length(idx)
        }

        # -- G observe par locus + permutations - moteur C++ natif (voir
        # src/g_test_subdivision.cpp, verifie ligne a ligne contre la logique
        # R ci-dessous), avec repli automatique sur R en cas d'echec.
        loci_idx0 <- lapply(loci_idx, function(ix) as.integer(ix - 1L))

        cpp_ok <- TRUE
        obs_res <- tryCatch(
          g_test_subdivision_observed_cpp(as.integer(pop_idx0_full), loci_idx0, loci_allele0,
                                           as.integer(loci_n_allele), n_pops),
          error = function(e) { cpp_ok <<- FALSE; NULL }
        )

        if (cpp_ok) {
          g_obs <- as.numeric(obs_res$g_obs); names(g_obs) <- loci_names
          g_obs_overall <- obs_res$g_obs_overall
        } else {
          # -- Repli R (identique a l'ancienne implementation) --------------
          g_obs <- vapply(seq_len(n_loci), function(j) {
            pop0_j <- pop_idx0_full[loci_idx[[j]]]
            if (loci_n_allele[j] < 2L || length(unique(pop0_j)) < 2L) return(NA_real_)
            .g_stat_vec(pop0_j, loci_allele0[[j]], n_pops, loci_n_allele[j])
          }, numeric(1))
          names(g_obs) <- loci_names
          g_obs_overall <- sum(g_obs, na.rm = TRUE)
        }


        # -- Permutations ------------------------------------------------------
        # H0 (FSTAT, NOT assuming HW within samples) : les GENOTYPES complets sont
        # reassignes au hasard entre populations. On tire une seule permutation
        # globale par replicat b (niveau individu), reutilisee pour tous les loci
        # (meme individu = meme re-affectation partout dans un meme replicat),
        # en ne gardant que les lignes valides de chaque locus.
        n_ind          <- length(pop_idx0_full)
        G_null_locus   <- matrix(NA_real_, nrow = n_perm, ncol = n_loci)
        G_null_overall <- numeric(n_perm)

        set.seed(as.integer(.seed()))

        # Moteur C++ : tourne par lots (jusqu'a 20) pour garder une barre de
        # progression fluide, tout en executant la boucle chaude (permutations
        # x loci) en code natif au lieu d'une boucle R + vapply().
        if (cpp_ok) {
          n_batches <- max(1L, min(20L, n_perm))
          base_size <- n_perm %/% n_batches; rem <- n_perm %% n_batches
          batch_sizes <- rep(base_size, n_batches)
          if (rem > 0L) batch_sizes[seq_len(rem)] <- batch_sizes[seq_len(rem)] + 1L
          batch_sizes <- batch_sizes[batch_sizes > 0L]

          done <- 0L
          cpp_ok <- tryCatch({
            for (bs in batch_sizes) {
              rb <- g_test_subdivision_batch_cpp(as.integer(pop_idx0_full), loci_idx0, loci_allele0,
                                                  as.integer(loci_n_allele), n_pops, bs)
              rows <- (done + 1L):(done + bs)
              G_null_locus[rows, ] <- rb$g_null_locus
              G_null_overall[rows] <- rb$g_null_overall
              done <- done + bs
            }
            TRUE
          }, error = function(e) FALSE)
        }

        # -- Repli R (identique a l'ancienne implementation) - utilise si le
        # moteur C++ (observe OU permutations) a echoue pour une raison
        # quelconque ; redemarre proprement avec la meme graine.
        if (!cpp_ok) {
          set.seed(as.integer(.seed()))
          for (b in seq_len(n_perm)) {
            perm_pop0 <- pop_idx0_full[sample.int(n_ind)]   # permutation globale des individus

            g_perm <- vapply(seq_len(n_loci), function(j) {
              if (loci_n_allele[j] < 2L) return(NA_real_)
              pop0_j <- perm_pop0[loci_idx[[j]]]
              if (length(unique(pop0_j)) < 2L) return(NA_real_)
              .g_stat_vec(pop0_j, loci_allele0[[j]], n_pops, loci_n_allele[j])
            }, numeric(1))

            G_null_locus[b, ]  <- g_perm
            G_null_overall[b]  <- sum(g_perm, na.rm = TRUE)

          }
        }


        # -- P-values : deux definitions, comme FSTAT (colonnes [>= obs] et [> obs]) --
        .pvals <- function(obs, null) {
          null <- null[is.finite(null)]
          if (!is.finite(obs) || length(null) == 0L) return(c(NA_real_, NA_real_))
          c(
            (pga_n_ge(null, obs) + 1) / (length(null) + 1),
            (pga_n_gt(null, obs) + 1) / (length(null) + 1)
          )
        }

        p_mat <- vapply(seq_len(n_loci), function(j) .pvals(g_obs[j], G_null_locus[, j]), numeric(2))
        p_ge_locus <- p_mat[1, ]; p_gt_locus <- p_mat[2, ]
        names(p_ge_locus) <- names(p_gt_locus) <- loci_names

        p_overall_vec <- .pvals(g_obs_overall, G_null_overall)
        p_ge_overall  <- p_overall_vec[1]
        p_gt_overall  <- p_overall_vec[2]

        # -- Tables finales ----------------------------------------------------
        # Ordre = ordre physique DuckDB (loci_names de hf_mat_r reordonne)
        per_locus_tbl <- data.frame(
          ID       = loci_names,
          N_geno   = loci_n_valid,
          G_obs    = g_obs,
          p_ge     = p_ge_locus,
          p_gt     = p_gt_locus,
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        overall_row <- data.frame(
          ID       = "Overall",
          N_geno   = NA_integer_,
          G_obs    = g_obs_overall,
          p_ge     = p_ge_overall,
          p_gt     = p_gt_overall,
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        # Overall toujours en derniere ligne - meme convention que FST
        final_tbl <- rbind(per_locus_tbl, overall_row)


        duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
        list(
          final_table    = final_tbl,
          g_obs_overall  = g_obs_overall,
          p_global       = p_ge_overall,
          p_global_gt    = p_gt_overall,
          duration       = duration,
          metadata       = list(
            n_perm       = n_perm,
            loci_names   = loci_names,
            pop_names    = pop_names,
            n_complete   = stats::setNames(as.integer(tabulate(pop_idx0_full + 1L, nbins = n_pops)), pop_names),
            dataset_name = if (!is.null(rv$dataset_filename)) rv$dataset_filename else NA_character_
          )
        )

      }, error = function(e) {
        showNotification(paste("Error in G-based test:", e$message), type = "error")
        NULL
      })

      ok
    }

    # ==================================== FIN G-TEST ===============================================

    ###

  })
}