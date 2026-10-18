# mod_global_panmixia.R

mod_global_panmixia_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("globe", "Global Panmixia · FIT",""),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    "FIT: CI & p-value parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm_fit"),    "Number of Permutations:",        value = 10000, min = 100, max = 20000, step = 100),
            numericInput(ns("n_boot_fit"),    "Number of Bootstrap Replicates:", value = 10000, min = 100, max = 20000, step = 100),
            numericInput(ns("conf_level_fit"),"Confidence Level:",               value = 0.95, min = 0.80, max = 0.99, step = 0.01),
            uiOutput(ns("ui_fit_out_status")),
            downloadButton(ns("Run_FIT_Analysis"), " Run",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;")
          ),
          column(9,
            h4(icon("chart-line"), "FIT Analysis Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3,
                valueBoxOutput(ns("global_fit_box"),       width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("global_fit_pvalue_box"),width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("significant_loci_fit_box"),  width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("analysis_time_fit_box"),width = NULL)
              )
            ),
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fit_progress"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    )
  )
}