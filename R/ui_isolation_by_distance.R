# ui_isolation_by_distance.R

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner("atom","Isolation by Distance · Mantel Test",""),

    fluidRow(
      box(
        width = 12, solidHeader = FALSE,
        uiOutput(ns("ui_run_status"))
      )
    ),

    # ══════════════════════════════════════════════════════════════════════
    # SECTION 1 — Isolation by Distance
    # ══════════════════════════════════════════════════════════════════════
    h2("Isolation by Distance", class = "section-title"),

    fluidRow(
      box(width = 12, solidHeader = TRUE, status = "primary",
          title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                      icon("database"), " Data source"),
        radioButtons(ns("ibd_source"), NULL,
          choices = c("Null Alleles module (this session)" = "internal",
                      "Re-load an exported pairwise file"   = "external"),
          selected = "internal", inline = TRUE),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'external'", ns("ibd_source")),
          fluidRow(
            column(6,
              fileInput(ns("ibd_ext_file"), "Pairwise file (Pop1, Pop2, FST_ENA, DCSE_INA, Dgeo_m\u2026):",
                        accept = c(".csv", ".txt", ".tsv")),
              radioButtons(ns("ibd_ext_sep"), "Separator:",
                choices = c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";"),
                selected = "\t", inline = TRUE),
              checkboxInput(ns("ibd_ext_header"), "File has header row", value = TRUE)
            ),
            column(6,
              uiOutput(ns("ibd_ext_file_status")),
              uiOutput(ns("ibd_ext_dgeo_warning"))
            )
          )
        )
      )
    ),

    fluidRow(
      box(width = 12, solidHeader = TRUE, status = "primary",
          title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                      icon("sliders-h"), " Rousset's Isolation by Distance model ",
                      tags$a(href = "#", style="font-weight:400;color:#4f7cff;font-size:0.85em;",
                        onclick = "var el=document.querySelector('a[data-value=\"help\"]'); if(el) el.click(); return false;",
                        "(please see documentation)")),
        fluidRow(
          column(6, uiOutput(ns("ibd_col_geo_ui"))),
          column(6, uiOutput(ns("ibd_col_avg_ui")))
        ),
        fluidRow(
          column(6,
            uiOutput(ns("ibd_col_lo_ui")),
            tags$p(style="color:#777;font-size:10.5px;margin-top:2px;", "Must be of the same format as the genetic distance used.")
          ),
          column(6,
            uiOutput(ns("ibd_col_hi_ui")),
            tags$p(style="color:#777;font-size:10.5px;margin-top:2px;", "Must be of the same format as the genetic distance used.")
          )
        ),
        tags$hr(),
        fluidRow(
          column(4,
            conditionalPanel(
              condition = sprintf("input['%s'] == 'external'", ns("ibd_source")),
              tags$p(style = "color:#777;font-size:11px;",
                "D_geo is read directly from the uploaded file's Dgeo_m / lnDgeo columns.")
            )
          ),
          column(4,
            tags$div(style="font-size:12px;color:#555;margin-bottom:4px;", "Output file name:"),
            textInput(ns("ibd_out_root"), NULL, value = "",
                      placeholder = "auto-filled from imported file"),
            tags$p(style="color:#777;font-size:11px;",
              "All output files will be saved in a zipped file.")
          ),
          column(4,
            tags$div(style="margin-top:22px;"),
            downloadButton(ns("run_ibd"), " Run",
                         icon = icon("rocket"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          )
        )
      )
    ),

    fluidRow(
      box(width = 12, solidHeader = TRUE, status = "primary",
          title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                      icon("chart-bar"), " Results"),
        uiOutput(ns("ui_ibd_status"))
      )
    ),

    # ══════════════════════════════════════════════════════════════════════
    # SECTION 2 — Mantel Test
    # ══════════════════════════════════════════════════════════════════════
    h2("Mantel Test", class = "section-title"),

    fluidRow(
      box(width = 12, solidHeader = TRUE, status = "primary",
          title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                      icon("database"), " Data source"),
        radioButtons(ns("mt_source"), NULL,
          choices = c("Internal pairwise table (from Null Alleles module)" = "internal",
                      "Upload a file"                                      = "upload"),
          selected = "internal", inline = TRUE),

        conditionalPanel(
          condition = sprintf("input['%s'] == 'upload'", ns("mt_source")),
          fluidRow(
            column(6,
              fileInput(ns("mt_file"), "File (Pop1, Pop2, dist1, dist2, ...):",
                        accept = c(".csv", ".txt", ".tsv")),
              radioButtons(ns("mt_sep"), "Separator:",
                choices = c("Tab"="\t", "Comma"=",", "Semicolon"=";"),
                selected = "\t", inline = TRUE),
              checkboxInput(ns("mt_header"), "File has header row", value = TRUE)
            ),
            column(6, uiOutput(ns("mt_file_status")))
          )
        ),

        conditionalPanel(
          condition = sprintf("input['%s'] == 'internal'", ns("mt_source")),
          checkboxInput(ns("mt_use_extra"),
            "Merge an extra distance file (temporal / ecological / categorical)",
            value = FALSE),
          conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("mt_use_extra")),
            fluidRow(
              column(6,
                fileInput(ns("mt_extra_file"), "Extra file (first 2 cols = Pop1, Pop2 IDs):",
                          accept = c(".csv", ".txt", ".tsv")),
                radioButtons(ns("mt_extra_sep"), "Separator:",
                  choices = c("Comma"=",", "Tab"="\t", "Semicolon"=";"),
                  selected = ",", inline = TRUE),
                checkboxInput(ns("mt_extra_header"), "File has header row", value = TRUE)
              ),
              column(6, uiOutput(ns("mt_extra_file_status")))
            )
          )
        ),

        tags$hr(),
        fluidRow(
          column(6, uiOutput(ns("mt_col_pop1_ui"))),
          column(6, uiOutput(ns("mt_col_pop2_ui")))
        )
      )
    ),

    fluidRow(
      box(width = 12, solidHeader = TRUE, status = "primary",
          title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                      icon("sliders-h"), " Mantel parameters"),
        fluidRow(
          column(3,
            checkboxGroupInput(ns("mt_stats"), "Statistic (select 1\u20134):",
              choices = c("Pearson r" = "r", "Spearman rho" = "spearman",
                          "Rousset's 1D" = "rousset1d", "Rousset's 2D" = "rousset2d"),
              selected = "r")
          ),
          column(3,
            radioButtons(ns("mt_p_formula"), "p-value formula:",
              choices = c("(b+1)/(m+1) \u2014 corrected" = "plus1",
                          "b/m \u2014 plain proportion"   = "plain"),
              selected = "plus1"),
            numericInput(ns("mt_n_perm"), "Permutations:",
                         value = 10000, min = 99, max = 200000, step = 1000)
          ),
          column(3,
            tags$div(style="font-size:12px;color:#555;margin-bottom:4px;", "Output file name:"),
            textInput(ns("mt_out_root"), NULL, value = "",
                      placeholder = "auto-filled from imported file"),
            tags$p(style="color:#777;font-size:11px;",
              "All output files will be saved in a zipped file.")
          ),
          column(3,
            downloadButton(ns("run_mantel"), " Run",
                         icon = icon("rocket"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          )
        ),
        tags$hr(),
        fluidRow(
          column(3,
            conditionalPanel(
              condition = sprintf("input['%s'] && (input['%s'].includes('r') || input['%s'].includes('spearman'))",
                                   ns("mt_stats"), ns("mt_stats"), ns("mt_stats")),
              tags$div(style="font-weight:600;color:#333a43;margin-bottom:6px;", "Pearson r / Spearman rho"),
              uiOutput(ns("mt_col_x_ui")),
              uiOutput(ns("mt_col_y_ui")),
              tags$p(style="color:#777;font-size:11px;", icon("lock"),
                " Same X and Y for both. A geographic distance column (raw or already-logged) is tested on the ",
                "log scale automatically; any other variable is used as-is \u2014 no setting needed.")
            )
          ),
          column(3,
            conditionalPanel(
              condition = sprintf("input['%s'] && input['%s'].includes('rousset1d')", ns("mt_stats"), ns("mt_stats")),
              tags$div(style="font-weight:600;color:#333a43;margin-bottom:6px;", "Rousset's 1D"),
              uiOutput(ns("mt_col_x_1d_ui")),
              uiOutput(ns("mt_col_y_1d_ui")),
              tags$p(style="color:#777;font-size:11px;", icon("lock"), " X is always the raw distance (D_geo).")
            )
          ),
          column(3,
            conditionalPanel(
              condition = sprintf("input['%s'] && input['%s'].includes('rousset2d')", ns("mt_stats"), ns("mt_stats")),
              tags$div(style="font-weight:600;color:#333a43;margin-bottom:6px;", "Rousset's 2D"),
              uiOutput(ns("mt_col_x_2d_ui")),
              uiOutput(ns("mt_col_y_2d_ui")),
              tags$p(style="color:#777;font-size:11px;", icon("lock"), " X is always ln(distance) \u2014 ln(D_geo).")
            )
          )
        )
      )
    ),

    fluidRow(
      box(width = 12, solidHeader = TRUE, status = "primary",
          title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                      icon("chart-bar"), " Results"),
        uiOutput(ns("ui_mantel_status"))
      )
    )
  )
}