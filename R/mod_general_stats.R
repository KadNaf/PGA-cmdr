# mod_general_stats.R

mod_general_stats_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    useWaiter(),
    tags$head(gs_head()),

    module_banner("table", "General Statistics",""),

    fluidRow(
      box(
        width = 6,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("chart-bar"), "Statistics selection"),
        solidHeader = TRUE, status = "primary",
        h4("Select Genetic Indices"),
        h5("Diversity Measures:"),
        checkboxInput(ns("ho_checkbox"),  "Ho (Observed Heterozygosity)", TRUE),
        checkboxInput(ns("hs_checkbox"),  "Hs (Expected Heterozygosity within populations)", TRUE),
        checkboxInput(ns("ht_checkbox"),  "Ht (Total expected heterozygosity)", TRUE),
        h5("F-statistics (Weir & Cockerham):"),
        checkboxInput(ns("fis_wc_checkbox"), "FIS (Weir & Cockerham estimator)", TRUE),
        checkboxInput(ns("fst_wc_checkbox"), "FST (Weir & Cockerham estimator)", TRUE),
        checkboxInput(ns("fit_wc_checkbox"), "FIT (Weir & Cockerham estimator)", TRUE),
        h5("Advanced Statistics:"),
        checkboxInput(ns("fst_max_checkbox"),  "Fst-max (Maximum differentiation, Meirmans)", FALSE),
        checkboxInput(ns("fst_prim_checkbox"), "Fst' (Meirmans) (Empirical standardisation)", FALSE),
        checkboxInput(ns("GST_checkbox"),      "GST (Nei's genetic differentiation)", FALSE),
        checkboxInput(ns("GST_sec_checkbox"),  "GST'' (Hedrick's correction)", FALSE),
        tags$hr(),
        uiOutput(ns("ui_gs_out_status")),
        downloadButton(ns("run_basic_stats"),
                     label = "Run",
                     icon = icon("rocket"),
                     class = "btn-action-primary btn-block")
      )
    )
  )
}