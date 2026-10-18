# server_import_data.R

server_import_data <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    rv$preview_raw  <- NULL
    rv$preview_meta <- NULL
    uploaded_file <- shiny::reactiveVal(NULL)
    
    
    
    # ---------------------------------------------------------#
    # Core formatting pipeline (shared by AUTO + MANUAL) ####
    # ---------------------------------------------------------#
    do_assign_and_format <- function(colnames_all,
                                     pop_data,
                                     col_ranges_data,
                                     metadata_ranges,
                                     missing_code,
                                     latitude_data,
                                     longitude_data,
                                     selected_levels = character(0),
                                     make_map = TRUE,
                                     con = NULL,
                                     tbl_raw = NULL) {
      
      # -------------------------#
      # Missing handling in duckDB
      # -------------------------#
      miss_info   <- normalize_missing_code(missing_code)
      missing_set <- miss_info$miss_set
      
      # -------------------------#
      # Required inputs
      # -------------------------#
      if (is.null(pop_data) || length(pop_data) != 1L || is.na(pop_data) || !nzchar(pop_data)) {
        shinyalert::shinyalert("Error", "You need to select populations.", type = "error")
        return(invisible(FALSE))
      }
      
      if (is.null(col_ranges_data) || length(col_ranges_data) != 1L || is.na(col_ranges_data) || !nzchar(col_ranges_data)) {
        shinyalert::shinyalert("Error", "You need to select a marker range.", type = "error")
        return(invisible(FALSE))
      }
      if (is.null(colnames_all) || !length(colnames_all)) {
        shinyalert::shinyalert("Error", "No column names available (header not read).", type = "error")
        return(invisible(FALSE))
      }
      
      n_all <- length(colnames_all)
      
      # -------------------------#
      # Marker range validation
      # -------------------------#
      rvv <- unlist(strsplit(col_ranges_data, "[:-]"))
      rvv <- suppressWarnings(as.integer(rvv))
      
      if (length(rvv) != 2 || any(is.na(rvv)) || rvv[1] < 1 || rvv[2] > n_all) {
        shinyalert::shinyalert("Error", "Try again, your range is out of bounds.", type = "error")
        return(invisible(FALSE))
      }
      
      marker_idx <- seq.int(min(rvv), max(rvv))
      marker_raw_names <- colnames_all[marker_idx]
      
      # -------------------------#
      # GPS detection
      # -------------------------#
      has_gps <- !is.null(latitude_data)  && nzchar(latitude_data) &&
        !is.null(longitude_data) && nzchar(longitude_data) &&
        latitude_data  %in% colnames_all &&
        longitude_data %in% colnames_all &&
        latitude_data != pop_data &&
        longitude_data != pop_data
      
      if (has_gps && (latitude_data %in% marker_raw_names || longitude_data %in% marker_raw_names)) {
        shinyalert::shinyalert("Error", "Latitude/Longitude cannot be inside the marker range.", type = "error")
        return(invisible(FALSE))
      }
      
      # -------------------------#
      # Metadata ranges
      # -------------------------#
      meta_idx   <- parse_col_index_ranges(metadata_ranges, n_max = n_all)
      meta_names <- if (length(meta_idx)) colnames_all[meta_idx] else character(0)
      
      # forbid overlaps
      meta_names <- setdiff(meta_names, c(pop_data, latitude_data, longitude_data, marker_raw_names))
      
      # -------------------------#
      # Build keep columns (names)
      # -------------------------#
      keep <- c(pop_data)
      if (has_gps) keep <- c(keep, latitude_data, longitude_data)
      if (length(meta_names)) keep <- c(keep, meta_names)
      
      # all markers in preview (rows limited to preview_n by SQL LIMIT)
      marker_preview_names <- marker_raw_names
      keep <- unique(c(keep, marker_preview_names))
      keep <- intersect(keep, colnames_all)
      
      if (!length(keep)) {
        shinyalert::shinyalert("Error", "No columns selected (check separator/header).", type = "error")
        return(invisible(FALSE))
      }

      # -------------------------#
      # Read only required columns
      # -------------------------#
      if (is.null(con) || is.null(tbl_raw) || !nzchar(tbl_raw)) {
        shinyalert::shinyalert("Error", "Database not initialised.", type = "error")
        return(invisible(FALSE))
      }
      
      ok_raw <- tryCatch(DBI::dbExistsTable(con, tbl_raw), error = function(e) FALSE)
      if (!isTRUE(ok_raw)) {
        shinyalert::shinyalert("Error", "Raw table not found in DuckDB (import step missing).", type = "error")
        return(invisible(FALSE))
      }
      
      preview_n <- 100L
      
      # quote identifiers safely (columns + table)
      keep_sql <- paste(
        vapply(keep, function(x) as.character(DBI::dbQuoteIdentifier(con, x)), character(1)),
        collapse = ", "
      )
      tbl_sql <- as.character(DBI::dbQuoteIdentifier(con, tbl_raw))
      
      # IMPORTANT: bring rowid as stable individual id (DuckDB)
      sql_full <- sprintf("SELECT rowid AS individual, %s FROM %s;", keep_sql, tbl_sql)
      
      new_df <- .time_it("dbGetQuery selected cols -> R", {
        tryCatch(DBI::dbGetQuery(con, sql_full), error = function(e) NULL)
      })
      if (is.null(new_df) || !nrow(new_df)) {
        shinyalert::shinyalert("Error", "No rows returned from DuckDB for selected columns.", type = "error")
        return(invisible(FALSE))
      }
      
      new_df <- as.data.frame(new_df, stringsAsFactors = FALSE)
      
      
      
      # -------------------------#
      # Standardize key column names
      # -------------------------#
      names(new_df)[names(new_df) == pop_data] <- "Population"
      
      if (has_gps) {
        names(new_df)[names(new_df) == latitude_data]  <- "Latitude"
        names(new_df)[names(new_df) == longitude_data] <- "Longitude"
        new_df$Latitude  <- suppressWarnings(as.numeric(new_df$Latitude))
        new_df$Longitude <- suppressWarnings(as.numeric(new_df$Longitude))
      }
      
      # -------------------------#
      # Marker frame
      # -------------------------#
      marker_raw_names_ok <- marker_raw_names
      if (!length(marker_raw_names_ok)) {
        shinyalert::shinyalert("Error", "No marker columns found (check marker range).", type = "error")
        return(invisible(FALSE))
      }
      
      # All markers in preview (rows already limited to 100 by LIMIT clause)
      marker_preview_names <- marker_raw_names_ok
      marker_preview_ok <- intersect(marker_preview_names, names(new_df))
      if (!length(marker_preview_ok)) {
        shinyalert::shinyalert("Error", "Marker preview columns not found in selected columns (check marker range).", type = "error")
        return(invisible(FALSE))
      }
      
      marker_df <- new_df[, marker_preview_ok, drop = FALSE]
      
      # -------------------------#
      # PREVIEW ONLY (NO DECODING)
      # -------------------------#
      
      meta_part <- .time_it("preview: extract meta part (R)", {
        new_df[, setdiff(names(new_df), marker_preview_ok), drop = FALSE]
      })
      
      marker_preview <- .time_it("preview: extract marker subset (R)", {
        new_df[, marker_preview_ok, drop = FALSE]
      })
      
      new_df_prev <- .time_it("preview: cbind meta + marker subset (R)", {
        out <- cbind(meta_part, marker_preview, stringsAsFactors = FALSE)
        
        # ---- merge paired allele columns into "a1/a2" for preview 
        # Example: B12 + B12_1  ->  B12 = "195/197"  (and remove B12_1)
        pick_b <- function(locus, nms) {
          cands <- c(paste0(locus, "_1"), paste0(locus, ".", 1:9))
          hit <- cands[cands %in% nms]
          if (length(hit)) hit[1] else NA_character_
        }
        
        nms <- names(out)
        loci <- unique(sub("(_1|\\.[0-9]+)$", "", nms))
        
        for (locus in loci) {
          a <- locus
          b <- pick_b(locus, nms)
          
          if (!a %in% nms) next
          if (is.na(b) || !b %in% nms) next
          
          out[[locus]] <- paste0(out[[a]], "/", out[[b]])
          
          # drop the allele columns (keep only merged locus)
          a_val <- as.character(out[[a]])
          b_val <- as.character(out[[b]])
          
          already_diploid <- grepl("/", a_val, fixed = TRUE) | grepl("-", a_val, fixed = TRUE)
          
          out[[a]] <- ifelse(
            already_diploid,
            a_val,                             
            paste0(a_val, "/", b_val)         
          )
          
          out[[b]] <- NULL 
          
          # refresh names after deletion
          nms <- names(out)
        }
        
        out
      })
      
      .time_it("preview: dbWriteTable formatted_preview", {
        DBI::dbWriteTable(con, rv$tbl_formatted_preview, new_df_prev, overwrite = TRUE)
      })
      
      base <- .time_it("base: infer from marker strings (R)", {
        infer_base_from_marker_strings(marker_df, sep = "/", default_base = 1000L)
      })
      haplotype_length <- suppressWarnings(as.integer(round(log10(as.numeric(base)))))
      if (!is.finite(haplotype_length) || haplotype_length < 1L) haplotype_length <- 3L
      
      .time_it("params: store downstream config", {
        .duckdb_set_params(con, list(
          is_preview       = 1L,
          preview_n        = as.integer(preview_n),
          pop_col_raw      = pop_data,
          missing_code_raw = as.character(missing_code),
          missing_gt       = 0L,
          missing_set      = paste(sort(unique(missing_set)), collapse = ","),
          sep              = "/",
          base             = as.integer(base),
          haplotype_length   = as.integer(haplotype_length),
          width_scalar_full  = as.integer(haplotype_length),
          base_scalar_full   = as.integer(base),
          marker_cols      = marker_raw_names_ok,   
          tbl_raw          = tbl_raw,
          tbl_meta         = rv$tbl_meta %||% "meta",
          tbl_formatted_preview = rv$tbl_formatted_preview %||% "formatted_preview",
          tbl_hf           = rv$tbl_hf %||% "hf",
          marker_cols_raw  = marker_raw_names_ok,                         # <-- CHANGE
          locus_cols       = .get_locus_cols_from_marker_cols(marker_raw_names_ok),
          genotype_format  = if (any(grepl("(_1|\\.[0-9]+)$", marker_raw_names_ok))) "paired" else "string"
        ))
      })
      
      .time_it("preview: formatted_preview from DB (LIMIT 100)", {
        rv$preview_raw <- DBI::dbGetQuery(
          con,
          sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_formatted_preview))
        )
      })
      
      rm(marker_df, meta_part, marker_preview, new_df_prev, new_df)
      gc(FALSE)
      
      invisible(TRUE)
    }
    
    
    # ---------------------------------------------------------#
    ## Session-scoped DuckDB ####
    # ---------------------------------------------------------#
    session_dir <- file.path(tempdir(), paste0("pgacmdr_", session$token))
    dir.create(session_dir, showWarnings = FALSE, recursive = TRUE)
    
    rv$tbl_raw               <- "raw"
    rv$tbl_meta              <- "meta"
    rv$tbl_formatted         <- "formatted"         # full (option B, later)
    rv$tbl_formatted_preview <- "formatted_preview" # UI-only
    rv$tbl_hf                <- "hf"
    rv$db_tick <- 0L
    
    session$userData$db_path <- file.path(session_dir, "pgacmdr.duckdb")
    
    con <- DBI::dbConnect(
      duckdb::duckdb(),
      dbdir = session$userData$db_path,
      read_only = FALSE
    )
    
    .duckdb_tune_for_big_import(con)
    
    rv$con <- con
    rv$db_ready <- TRUE

    
    
    session$onSessionEnded(function() {
      try(DBI::dbDisconnect(rv$con, shutdown = TRUE), silent = TRUE)
      try(duckdb::duckdb_shutdown(), silent = TRUE)
      try(unlink(session$userData$db_path, force = TRUE), silent = TRUE)
      if (dir.exists(session_dir)) unlink(session_dir, recursive = TRUE, force = TRUE)
    })
    
    # ---------------------------------------------------------#
    # Load_default_data ####
    # ---------------------------------------------------------#
    
    shiny::observeEvent(input$load_default_data, {
      
      reset_downstream_state(rv)
      
      default_path <- system.file("extdata", "default_dataset.csv", package = "pgacmdr")
      if (!nzchar(default_path) || !file.exists(default_path)) {
        shinyalert::shinyalert("Error", "Default dataset not found in package (inst/extdata/default_dataset.csv).", type = "error")
        return()
      }

      rv$file_path <- default_path
      rv$sep <- "\t"
      rv$header <- input$header
      rv$dataset_filename <- "default_dataset.csv"
      
      # ------------------------------#
      # State 1: import RAW
      # ------------------------------#
      ok_db <- tryCatch(
        .duckdb_import_raw(
          con       = rv$con,
          tbl_raw   = rv$tbl_raw,
          file_path = rv$file_path,
          sep       = rv$sep,
          header    = rv$header
        ),
        error = function(e) {
          shinyalert::shinyalert("DuckDB import failed", conditionMessage(e), type = "error")
          FALSE
        }
      )
      if (!isTRUE(ok_db)) return()
      
      # ------------------------------#
      # State 2: header + raw preview from DuckDB
      # ------------------------------#
      rv$colnames_all <- tryCatch(DBI::dbListFields(rv$con, rv$tbl_raw), error = function(e) character(0))
      if (!length(rv$colnames_all)) {
        shinyalert::shinyalert("Header read failed", "No columns detected in DuckDB raw table.", type = "error")
        return()
      }
      
      update_metadata_choices(session, rv$colnames_all)
      
      rv$preview_raw <- tryCatch(
        DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_raw))),
        error = function(e) NULL
      )
      if (is.null(rv$preview_raw) || !nrow(rv$preview_raw)) {
        shinyalert::shinyalert("No data rows", "Raw table has 0 rows after import.", type = "warning")
        return()
      }
      
      # ------------------------------#
      # State 3: auto-detect on raw preview
      # ------------------------------#
      miss_auto <- normalize_missing_code(input$missing_code)
      det <- detect_columns_auto(rv$preview_raw, missing_info = miss_auto)
      rv$det <- det
      
      # populate the manual panel selections
      populate_manual_from_detection(session, det, rv$colnames_all)
      
      # meta_ranges used BOTH for meta table + for auto-format call
      meta_idx <- match(det$metadata_cols, rv$colnames_all)
      meta_idx <- meta_idx[!is.na(meta_idx)]
      meta_ranges <- .compress_idx_ranges(meta_idx)
      
      # If population missing -> stop (manual UI is ready)
      if (is.null(det$population) || is.na(det$population) || !nzchar(det$population)) {
        shinyalert::shinyalert(
          "Auto-detect incomplete",
          "Population column not detected. Either change the separator or use the manual mode.",
          type = "warning"
        )
        return()
      }
      
      # ------------------------------#
      # State 3b: AUTO-FORMAT 
      # ------------------------------#
      if (!is.null(det$population) && nzchar(det$population) &&
          !is.null(det$marker_range) && nzchar(det$marker_range)) {
        
        ok_fmt <- do_assign_and_format(
          colnames_all    = rv$colnames_all,
          pop_data        = det$population,
          col_ranges_data = det$marker_range,
          metadata_ranges = meta_ranges,
          missing_code    = input$missing_code,
          latitude_data   = det$latitude %|||% "",
          longitude_data  = det$longitude %|||% "",
          selected_levels = character(0),
          make_map        = FALSE,
          con             = rv$con,
          tbl_raw         = rv$tbl_raw
        )
        
        if (isTRUE(ok_fmt)) {
          
          # full DB: formatted is already written by do_assign_and_format()
          ok_fmt_tbl <- tryCatch(DBI::dbExistsTable(rv$con, rv$tbl_formatted_preview), error = function(e) FALSE)
          if (!isTRUE(ok_fmt_tbl)) {
            shinyalert::shinyalert("Error", "Formatted table not found in DuckDB.", type = "error")
            return()
          }
        }
      }
      
      # ------------------------------#
      # State 4: build META + preview + map (DB meta used for map)
      # ------------------------------#
      tryCatch(
        .time_it("META build", {
          .duckdb_build_meta(
            con = rv$con,
            tbl_raw = rv$tbl_raw,
            tbl_meta = rv$tbl_meta,
            colnames_all = rv$colnames_all,
            pop_data = det$population,
            latitude_data = det$latitude,
            longitude_data = det$longitude,
            metadata_ranges = meta_ranges
          )
        })
        ,
        error = function(e) {
          shinyalert::shinyalert("DuckDB meta build failed", conditionMessage(e), type = "error")
          NULL
        }
      )
      
      # ---- ensure meta.individual exists and matches raw.rowid 
      DBI::dbExecute(
        rv$con,
        sprintf(
          "
    ALTER TABLE %s ADD COLUMN IF NOT EXISTS individual BIGINT;

    UPDATE %s m
    SET individual = r.rowid
    FROM %s r
    WHERE m.rowid = r.rowid;
    ",
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_raw)
        )
      )
      
    # ---- build HF (chunked) AFTER meta exists
      .build_hf_from_params <- function(con, rv, missing_code_raw, batch_size = 10000L) {
        shiny::req(con)
        
        # marker_cols_raw (raw columns, may include suffixes like _1 / .1)
        marker_json <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
          error = function(e) NA_character_
        )
        
        marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json)) {
          tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
        } else character(0)
        
        if (!length(marker_cols_raw)) stop("No marker_cols_raw found in DuckDB params.")
        
        # base
        base <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='base'")$value[1],
          error = function(e) NA_character_
        )
        base <- suppressWarnings(as.integer(base))
        if (!is.finite(base) || base <= 0L) stop("Invalid base in DuckDB params.")
        
        # genotype format
        fmt <- .duckdb_get_param(con, "genotype_format", default = "auto")
        if (identical(fmt, "auto")) {
          fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"
        }
        
        .time_it("HF build (chunked)", {
          .duckdb_build_hf_from_raw_chunked(
            con              = con,
            tbl_raw          = rv$tbl_raw,
            tbl_meta         = rv$tbl_meta,
            tbl_hf           = rv$tbl_hf,
            marker_cols      = marker_cols_raw,   # <-- FIX
            genotype_format  = fmt,
            base             = base,
            missing_code_raw = missing_code_raw,
            missing_gt       = 0L,
            batch_size       = as.integer(batch_size),
            attach_pop_code  = TRUE
          )
        })
        
        invisible(TRUE)
      }
      
      tryCatch(
        {
          .build_hf_from_params(rv$con, rv, missing_code_raw = input$missing_code, batch_size = 10000L)
          
          rv$preview_meta <- tryCatch(
            DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_meta))),
            error = function(e) NULL
          )
          
          rv$db_tick <- rv$db_tick + 1L
        },
        error = function(e) {
          shinyalert::shinyalert("HF build failed", conditionMessage(e), type = "error")
        }
      )
    
      
    })
    
    
    # ---------------------------------------------------------#
    # Upload user data ####
    # ---------------------------------------------------------#
    shiny::observeEvent(input$file1, {
      shiny::req(input$file1)
      shiny::req(input$file1$datapath)

      uploaded_file(input$file1$datapath[1])

      # Keep the ORIGINAL file name (not the temp datapath) so that other
      # modules (e.g. Subdivision, Isolation by distance) can stamp exported
      # result files with the name of the source dataset.
      rv$dataset_filename <- input$file1$name[1]

      reset_downstream_state(rv)

      rv$file_path <- uploaded_file()
      rv$sep <- input$sep
      rv$header <- input$header
      # ---------------------------------------------------------#
      # State 1: import RAW
      # ---------------------------------------------------------#

      ok_db <- tryCatch(
        .duckdb_import_raw(
          con       = rv$con,
          tbl_raw   = rv$tbl_raw,
          file_path = rv$file_path,
          sep       = rv$sep,
          header    = rv$header
        ),
        error = function(e) {
          shinyalert::shinyalert("DuckDB import failed", conditionMessage(e), type = "error")
          FALSE
        }
      )
      if (!isTRUE(ok_db)) return()
      # ---------------------------------------------------------#
      # State 2: header + raw preview
      # ---------------------------------------------------------#

      rv$colnames_all <- tryCatch(DBI::dbListFields(rv$con, rv$tbl_raw), error = function(e) character(0))
      if (!length(rv$colnames_all)) {
        shinyalert::shinyalert("Header read failed", "No columns detected in DuckDB raw table.", type = "error")
        return()
      }

      # No automatic column recognition: only populate the dropdown choices.
      # The operator picks Population / Latitude / Longitude and the loci
      # range themselves, then clicks "Apply" (see run_assign below).
      populate_choices_only(session, rv$colnames_all)

      rv$preview_raw <- tryCatch(
        DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_raw))),
        error = function(e) NULL
      )
      if (is.null(rv$preview_raw) || !nrow(rv$preview_raw)) {
        shinyalert::shinyalert("No data rows", "Raw table has 0 rows after import.", type = "warning")
        return()
      }

      rv$load_status <- sprintf(
        "File loaded: %d rows \u00d7 %d columns. Choose Population, Latitude/Longitude (optional) and the loci range below, then click \u201cApply\u201d.",
        nrow(rv$preview_raw), length(rv$colnames_all)
      )
    }, ignoreInit = TRUE)

    output$ui_load_status <- renderUI({
      shiny::req(rv$load_status)
      tags$p(style = "color:#166534;background:#f0fdf4;border:1px solid #bbf7d0;border-radius:4px;padding:6px 8px;font-size:12px;margin-top:8px;",
        icon("check-circle"), " ", rv$load_status)
    })

    # ── Numbered column reference: lets the operator see exactly which
    #    number corresponds to which column name, so the loci range (and
    #    Population/Latitude/Longitude choices) can be entered correctly.
    output$ui_columns_index <- renderUI({
      cn <- rv$colnames_all
      shiny::req(length(cn) > 0)
      df <- data.frame(Index = seq_along(cn), Column = cn, check.names = FALSE)
      n <- nrow(df)
      half <- ceiling(n / 2)
      left  <- df[seq_len(half), , drop = FALSE]
      right <- if (n > half) df[(half + 1):n, , drop = FALSE] else df[0, , drop = FALSE]

      .col_table <- function(d) {
        if (!nrow(d)) return(NULL)
        tags$table(class = "table table-striped", style = "width:100%; font-size:12px;",
          tags$thead(tags$tr(
            tags$th(style = "text-align:right; width:30%;", "Index"),
            tags$th(style = "text-align:left;", "Column")
          )),
          tags$tbody(
            lapply(seq_len(nrow(d)), function(i) {
              tags$tr(
                tags$td(style = "text-align:right;", d$Index[i]),
                tags$td(style = "text-align:left;", d$Column[i])
              )
            })
          )
        )
      }

      fluidRow(
        column(6, .col_table(left)),
        column(6, .col_table(right))
      )
    })

    # ── Loci range: Create-style input (number of loci + first locus
    #    column), with a read-only preview of the resulting range kept for
    #    verification. ──────────────────────────────────────────────────────
    # ── Detect whether loci are stored as ONE combined genotype column
    #    each, or as a PAIR of columns (allele1 in the base name, allele2
    #    in a "<name>_1" / "<name>.1" column right after it) — very common
    #    in raw genotype files. If the data is paired, "N loci" must span
    #    2*N raw columns, not N.
    .cols_per_locus <- function(colnames_all, first_col) {
      n_all <- length(colnames_all)
      if (is.na(first_col) || first_col < 1 || first_col >= n_all) return(1L)
      name1 <- colnames_all[first_col]
      name2 <- colnames_all[first_col + 1L]
      if (identical(name2, paste0(name1, "_1")) || identical(name2, paste0(name1, ".1"))) 2L else 1L
    }

    # ── Auto-suggest "Number of loci": as soon as the operator enters the
    #    first locus column, compute how many complete loci fit from there
    #    to the end of the file (accounting for paired-column format), and
    #    pre-fill the field with that suggestion — the operator can still
    #    adjust it down if some of those trailing columns aren't loci.
    observeEvent(input$first_locus_col, {
      first_col <- suppressWarnings(as.integer(input$first_locus_col))
      cn <- rv$colnames_all
      n_all <- length(cn %||% character(0))
      if (is.na(first_col) || first_col < 1 || first_col > n_all || n_all == 0) return(invisible(NULL))

      cpl <- .cols_per_locus(cn, first_col)
      n_suggested <- (n_all - first_col + 1L) %/% cpl
      if (n_suggested >= 1L) {
        updateNumericInput(session, "n_loci", value = n_suggested)
      }
    }, ignoreInit = TRUE)

    output$ui_locus_range_preview <- renderUI({
      n_loci <- suppressWarnings(as.integer(input$n_loci))
      first_col <- suppressWarnings(as.integer(input$first_locus_col))
      n_all <- length(rv$colnames_all %||% character(0))

      if (is.na(n_loci) || is.na(first_col) || n_loci < 1 || first_col < 1) {
        return(tags$p(style = "color:#777;font-size:13px;",
          "Enter the number of loci and the column number of the first locus."))
      }

      cpl <- .cols_per_locus(rv$colnames_all, first_col)
      last_col <- first_col + (n_loci * cpl) - 1L
      out_of_bounds <- n_all > 0 && (first_col > n_all || last_col > n_all)
      math_line <- tags$p(style = "margin:4px 0 0 0; color:#555;",
        sprintf("%d loci \u00d7 %d column%s/locus = columns %d to %d.",
                n_loci, cpl, if (cpl > 1) "s" else "", first_col, last_col))

      if (out_of_bounds) {
        tags$div(style = "color:#92400e;background:#fffbeb;border:1px solid #fcd34d;border-radius:4px;padding:8px 10px;font-size:13px;",
          tags$p(style = "margin:0;",
            icon("exclamation-triangle"),
            sprintf(" Out of bounds: this needs columns %d\u2013%d, but the file only has %d columns. Please try again.",
                    first_col, last_col, n_all)),
          math_line)
      } else {
        tags$div(style = "color:#166534;background:#f0fdf4;border:1px solid #bbf7d0;border-radius:4px;padding:8px 10px;font-size:13px;",
          tags$p(style = "margin:0;",
            icon("check-circle"), " Verification: ", tags$strong(sprintf("%d loci", n_loci)),
            if (cpl == 2L) " (paired-column format detected)" else ""),
          math_line)
      }
    })

    # ---------------------------------------------------------#
    shiny::observeEvent(input$run_assign, {
      shiny::req(rv$file_path, rv$colnames_all)

      # Loci range: computed from "number of loci" + "first locus column"
      # (Create-style), re-validated here exactly like the live preview.
      n_loci    <- suppressWarnings(as.integer(input$n_loci))
      first_col <- suppressWarnings(as.integer(input$first_locus_col))
      n_all     <- length(rv$colnames_all)

      if (is.na(n_loci) || is.na(first_col) || n_loci < 1 || first_col < 1) {
        shinyalert::shinyalert("Error", "Enter the number of loci and the column number of the first locus.", type = "error")
        return()
      }
      cpl_check <- .cols_per_locus(rv$colnames_all, first_col)
      last_col <- first_col + (n_loci * cpl_check) - 1L
      if (first_col > n_all || last_col > n_all) {
        msg <- sprintf("This needs columns %d\u2013%d (%d loci \u00d7 %d column%s/locus), but the file only has %d columns.",
                        first_col, last_col, n_loci, cpl_check, if (cpl_check > 1) "s" else "", n_all)
        shinyalert::shinyalert("Error", msg, type = "error")
        return()
      }
      col_ranges_data <- sprintf("%d-%d", first_col, last_col)

      # No "metadata columns" concept anymore: only Population / Latitude /
      # Longitude / loci are kept.
      meta_ranges <- ""

      ok <- do_assign_and_format(
        colnames_all    = rv$colnames_all,
        pop_data        = input$pop_data,
        col_ranges_data = col_ranges_data,
        metadata_ranges = meta_ranges,
        missing_code    = input$missing_code,
        latitude_data   = input$latitude_data,
        longitude_data  = input$longitude_data,
        make_map        = FALSE,
        con             = rv$con,
        tbl_raw         = rv$tbl_raw
      )
      if (!isTRUE(ok)) return()
      
      rv$preview_raw <- DBI::dbGetQuery(
        rv$con,
        sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_formatted_preview))
      )
      
      # ------------------------------#
      # State: build META 
      # ------------------------------#
      
      tryCatch(
        .time_it("META build", {
          .duckdb_build_meta(
            con             = rv$con,
            tbl_raw         = rv$tbl_raw,
            tbl_meta        = rv$tbl_meta,
            colnames_all    = rv$colnames_all,
            pop_data        = input$pop_data,
            latitude_data   = input$latitude_data,
            longitude_data  = input$longitude_data,
            metadata_ranges = meta_ranges
          )
        }),
        error = function(e) {
          shinyalert::shinyalert("DuckDB meta build failed", conditionMessage(e), type = "error")
          NULL
        }
      )
      # ---- ensure meta.individual exists and matches raw.rowid
      DBI::dbExecute(
        rv$con,
        sprintf(
          "
  ALTER TABLE %s ADD COLUMN IF NOT EXISTS individual BIGINT;

  UPDATE %s m
  SET individual = r.rowid
  FROM %s r
  WHERE m.rowid = r.rowid;
",
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_raw)
        )
      )
      
      # ---- build HF (chunked) AFTER meta exists
      .build_hf_from_params <- function(con, rv, missing_code_raw, batch_size = 10000L) {
        shiny::req(con)
        
        # marker_cols_raw (raw columns, may include suffixes like _1 / .1)
        marker_json <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
          error = function(e) NA_character_
        )
        
        marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json)) {
          tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
        } else character(0)
        
        if (!length(marker_cols_raw)) stop("No marker_cols_raw found in DuckDB params.")
        
        # base
        base <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='base'")$value[1],
          error = function(e) NA_character_
        )
        base <- suppressWarnings(as.integer(base))
        if (!is.finite(base) || base <= 0L) stop("Invalid base in DuckDB params.")
        
        # genotype format
        fmt <- .duckdb_get_param(con, "genotype_format", default = "auto")
        if (identical(fmt, "auto")) {
          fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"
        }
        
        .time_it("HF build (chunked)", {
          .duckdb_build_hf_from_raw_chunked(
            con              = con,
            tbl_raw          = rv$tbl_raw,
            tbl_meta         = rv$tbl_meta,
            tbl_hf           = rv$tbl_hf,
            marker_cols      = marker_cols_raw,   # <-- FIX
            genotype_format  = fmt,
            base             = base,
            missing_code_raw = missing_code_raw,
            missing_gt       = 0L,
            batch_size       = as.integer(batch_size),
            attach_pop_code  = TRUE
          )
        })
        
        invisible(TRUE)
      }
      
      tryCatch(
        {
          .build_hf_from_params(rv$con, rv, missing_code_raw = input$missing_code, batch_size = 10000L)
          
          rv$preview_meta <- tryCatch(
            DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_meta))),
            error = function(e) NULL
          )
          
          rv$db_tick <- rv$db_tick + 1L
          rv$load_status <- "Computations completed."
        },
        error = function(e) {
          shinyalert::shinyalert("HF build failed", conditionMessage(e), type = "error")
        }
      )
    })
  })
}
