# mod_local_panmixia.R

mod_local_panmixia_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("flask", "Local Panmixia \u00b7 FIS",""),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    "FIS: CI & p-value parameters (one Run = one file with 4 tables)"),
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
                         class = "btn-action-primary btn-block", style = "font-weight: bold;"),
            run_folder_warning()
          ),
          column(9,
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
    )
  )
}