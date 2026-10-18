# mod_genetic_diversities.R

mod_genetic_diversities_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("chart-line", "Genetic Diversities · HS · HT",""),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    "Genetic Diversity Analysis parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm_fst_div"),    "Number of Permutations:",        value = 10000, min = 100, max = 20000, step = 100),
            numericInput(ns("n_boot_fst_div"),    "Number of Bootstrap Replicates:", value = 10000, min = 100, max = 20000, step = 100),
            checkboxInput(ns("div_detail"), "Add the bootstrap replicate values (detail) at the end of the file", value = TRUE),
            uiOutput(ns("ui_div_out_status")),
            downloadButton(ns("run_FST_Analysis_div"), " Run",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;")
          ),
          column(9,
            h4(icon("chart-line"), "Diversity Analysis Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3, valueBoxOutput(ns("global_fst_div_box"),    width = NULL)),
              column(3, valueBoxOutput(ns("global_hs_box"),         width = NULL)),
              column(3, valueBoxOutput(ns("global_ht_box"),         width = NULL)),
              column(3, valueBoxOutput(ns("analysis_time_div_box"), width = NULL))
            ),
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fst_progress_div"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    )
  )
}