# mod_local_panmixia.R

mod_local_panmixia_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("flask", "Local Panmixia · FIS",""),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    "FIS: CI & p-value parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm"),    "Number of Permutations:",        value = 10000, min = 100, max = 20000),
            numericInput(ns("n_boot"),    "Number of Bootstrap Replicates:", value = 10000, min = 100, max = 20000),
            numericInput(ns("conf_level"),"Confidence Level:",               value = 0.95, min = 0.80, max = 0.99, step = 0.01),
            uiOutput(ns("ui_fis_out_status")),
            downloadButton(ns("Run_FIS_Analysis"), " Run",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;")
          ),
          column(9,
            h4(icon("chart-line"), "Analysis Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3, valueBoxOutput(ns("global_fis_box"),        width = NULL)),
              column(3, valueBoxOutput(ns("global_pvalue_box"),     width = NULL)),
              column(3, valueBoxOutput(ns("significant_loci_box"),  width = NULL)),
              column(3, valueBoxOutput(ns("analysis_time_box"),     width = NULL))
            ),
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fis_progress"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    ),

    h2("FIS · By Locus \u00d7 Population", class = "section-title"),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("table"),
                    "FIS (WC84) per locus \u00d7 population"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("fis_lp_n_perm"), "Number of Permutations:",
                         value = 10000, min = 100, max = 20000)
          ),
          column(9,
            br(),
            uiOutput(ns("ui_fislp_out_status")),
            downloadButton(ns("run_fis_locus_pop"), " Run",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;")
          )
        ),
        style = "padding: 10px;"
      )
    )
  )
}