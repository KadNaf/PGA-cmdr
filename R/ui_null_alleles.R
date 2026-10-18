# module/ui_null_alleles.R

null_alleles_UI <- function(id) {
  ns <- NS(id)

  # ── Supplemental CSS — only the bits with no shared-system equivalent:
  #    per-locus coding grid, pairwise matrices, bootstrap-result callouts.
  #    Recolored to the app's own palette (no more standalone dark/neon theme).
  supplemental_css <- tags$style(HTML("
    .na-info {
      background:#f5f7fa; border-left:4px solid #8D8680; border-radius:3px;
      padding:10px 14px; font-size:13px; line-height:1.65; color:#2c3e50; margin-bottom:14px;
    }
    .na-warn {
      background:#fff8e6; border-left:4px solid #E1AF00; border-radius:3px;
      padding:10px 14px; font-size:13px; line-height:1.65; color:#5c4400; margin-bottom:14px;
    }
    .na-locus-grid { display:flex; flex-wrap:wrap; gap:8px; margin-top:8px; }
    .na-locus-item {
      background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px;
      padding:.5rem .8rem; min-width:150px; flex:1;
    }
    .na-locus-item .control-label { display:none; } /* hide redundant label */
    .na-locus-name {
      font-size:20px; font-weight:700; color:#333a43;
      font-family:'Consolas','IBM Plex Mono',monospace; margin-bottom:4px;
    }
    .na-locus-item .radio { margin:2px 0; }
    .na-locus-item .radio label { font-size:13px; color:#333a43; }

    .na-boot-result {
      background:#f5f7fa; border:1px solid #8ea1b9; border-radius:6px;
      padding:.65rem 1rem; font-size:12.5px; color:#333a43;
      font-family:'Consolas','IBM Plex Mono',monospace; line-height:1.9; margin-top:.5rem;
    }
    .na-boot-result strong { color:#0c4a6e; }

    .na-matrix-wrap { overflow-x:auto; margin-top:.5rem; }
    .na-matrix { border-collapse:collapse; font-size:11px; font-family:'Consolas','IBM Plex Mono',monospace; width:100%; }
    .na-matrix th { background:#f8fafc; color:#475569; font-weight:600; padding:4px 9px; border:1px solid #e2e8f0; font-size:10.5px; white-space:nowrap; }
    .na-matrix td { padding:4px 9px; border:1px solid #e2e8f0; color:#1e293b; text-align:right; white-space:nowrap; font-size:11px; }
    .na-matrix tr:nth-child(even) td { background:#f8fafc; }
    .na-matrix .diag  { background:#f1f5f9 !important; color:#94a3b8; text-align:center; }
    .na-matrix .upper { color:#cbd5e1; text-align:center; }
    .na-matrix .lbl   { font-weight:700; color:#333a43; text-align:left; white-space:nowrap; }

    .na-dl-row { display:flex; gap:6px; flex-wrap:wrap; margin-top:.5rem; }
    .na-dl-row .btn { font-size:11px; padding:3px 12px; }

    .na-filecard { margin-bottom: 14px; }
    .na-filecard .fname { margin-top:8px; font-size:11px; word-break:break-all; color:#666; }
  "))

  fluidPage(
    tags$head(gs_head()),
    supplemental_css,

    module_banner("circle-notch", "Null Allele Estimation · FST-ENA · DCSE-INA",""),

    # ════════════════════════════════════════════════════════════════════
    # SETUP
    # ════════════════════════════════════════════════════════════════════
    fluidRow(
      box(
        width = 12,
        solidHeader = TRUE, status = "primary",
        tags$div(class = "na-warn",
          tags$p(style = "margin:.25rem 0;",
            "Choose the missing data code for each locus: ",
            tags$strong("0"), " = true missing (ignored); ",
            tags$strong("999999"), " = homozygote for allele 999 (null alleles)."
          ),
          tags$p(
            style = "margin:.25rem 0 0;font-weight:600;",
            "Ensure allele 999 is not already present in your dataset."
          )
        ),

        uiOutput(ns("locus_coding_ui")),

        tags$hr(),

        h4(icon("dice"), "Bootstrap parameters"),
        p("Bootstrap over subsamples (for 5 subsamples at least) and over loci (for 5 loci at least)"),
        tags$div(style = "display:flex; align-items:flex-start; gap:0; flex-wrap:wrap;",
          tags$div(style = "flex:1; min-width:190px; padding-right:16px;",
            numericInput(ns("nboot"),
              label = "Bootstraps over loci (at least 100):",
              value = 5000, min = 100, max = 99999, step = 1000, width = "100%")),
          tags$div(style = "border-left:1px solid #dcdfe4; flex:1; min-width:190px; padding:0 16px;",
            numericInput(ns("nboot_subs"),
              label = "Bootstraps over subsamples (at least 100):",
              value = 5000, min = 100, max = 99999, step = 1000, width = "100%"))
        ),
        tags$div(style = "display:flex; align-items:flex-start; gap:0; flex-wrap:wrap; margin-top:8px;",
          tags$div(style = "flex:1; min-width:190px; padding-right:16px;",
            numericInput(ns("alpha"),
              label = "Alpha (between 0.9999 and 0.0001):",
              value = 0.05, min = 0.0001, max = 0.5, step = 0.01, width = "100%")),
          tags$div(style = "border-left:1px solid #dcdfe4; flex:1; min-width:190px; padding:0 16px;",
            numericInput(ns("boot_seed"),
              label = "Random seed (between 1 and 1000000):",
              value = 12345, min = 1, max = 1000000, step = 1, width = "100%"))
        ),

        tags$hr(),

        h4(icon("map-marker-alt"), "GPS coordinates of subsamples"),
        radioButtons(ns("gps_available"), "Do you have GPS coordinates for your subsamples (decimal degrees)?",
          choices = c("Yes" = "yes", "No" = "no"), selected = "yes", inline = TRUE),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'yes'", ns("gps_available")),
          fluidRow(
            column(3, uiOutput(ns("gps_lon_col_ui"))),
            column(3, uiOutput(ns("gps_lat_col_ui")))
          )
        ),

        tags$hr(),

        h4(icon("save"), "Choose a name for the output"),
        tags$div(style = "max-width:320px;",
          textInput(ns("out_root"), NULL, value = "", placeholder = "auto-filled from the imported data file name")),
        tags$p(style="color:#777;font-size:11px;",
          "All output files will be saved in a zipped file."),

        tags$hr(),

        fluidRow(
          column(4,
            downloadButton(ns("run_all"),
              label = " Run",
              icon = icon("rocket"),
              class = "btn-action-primary btn-block",
              style = "font-weight: bold;"))
        ),
        br(),
        uiOutput(ns("ui_run_status"))
      )
    ),

    # ════════════════════════════════════════════════════════════════════
    # OUTPUT FILES — one card per exported file
    # ════════════════════════════════════════════════════════════════════
    h2("Output files", class = "section-title"),
    fluidRow(
      box(
        width = 12,
        title = uiOutput(ns("ui_output_files_title"), inline = TRUE),
        solidHeader = TRUE, status = "primary",
        tags$div(class = "pga-module-card na-filecard", style = "margin-bottom:14px; max-width:400px;",
          h5("p_nulls / locus"),
          p("Null allele frequencies per locus and subsamples, and averaged over subsamples."),
          tags$div(class = "fname", uiOutput(ns("ui_filename_1"), inline = TRUE)),
          uiOutput(ns("ui_dl_file1"))
        ),
        tags$div(class = "pga-module-card na-filecard", style = "margin-bottom:14px; max-width:400px;",
          h5("FST / FST-ENA"),
          p("Global FST per locus and over all, corrected or not for null alleles, with CI of bootstrap over subsamples and loci."),
          tags$div(class = "fname", uiOutput(ns("ui_filename_2"), inline = TRUE)),
          uiOutput(ns("ui_dl_file2"))
        ),
        tags$div(class = "pga-module-card na-filecard", style = "margin-bottom:14px; max-width:400px;",
          h5("Per-locus half-matrices"),
          p("Paired genetic distances in half left matrices, for use by other software."),
          tags$div(class = "fname", uiOutput(ns("ui_filename_4"), inline = TRUE)),
          uiOutput(ns("ui_dl_file4"))
        ),
        tags$div(class = "pga-module-card na-filecard", style = "margin-bottom:14px; max-width:400px;",
          h5("Bootstrap distributions"),
          p("Detailed bootstrap distribution for global FST's."),
          tags$div(class = "fname", uiOutput(ns("ui_filename_5"), inline = TRUE)),
          uiOutput(ns("ui_dl_file5"))
        ),
        tags$div(class = "pga-module-card na-filecard", style = "margin-bottom:14px; max-width:400px;",
          h5("Run parameters"),
          p("List of parameters you chose to use."),
          tags$div(class = "fname", uiOutput(ns("ui_filename_6"), inline = TRUE)),
          uiOutput(ns("ui_dl_file6"))
        ),
        uiOutput(ns("ui_file7_card"))
      )
    )
  )
}