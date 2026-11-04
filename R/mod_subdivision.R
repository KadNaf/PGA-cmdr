# mod_subdivision_ui.R

mod_subdivision_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("sitemap", "Population Subdivision, General Stats and standard errors",""),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    "FST and G-based test: CI & p-value parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm_fst"),     "Number of Permutations (FST and G-test):",        value = 10000, min = 100,  max = 20000, step = 100),
            numericInput(ns("n_boot_fst"),     "Number of Bootstrap Replicates:", value = 10000, min = 100,  max = 20000, step = 100),
            numericInput(ns("conf_level_fst"), "Confidence Level:",               value = 0.95, min = 0.80, max = 0.99,  step = 0.01),
            checkboxInput(ns("fst_detail"), "Add the bootstrap replicate values (detail) at the end of the file", value = TRUE),
            uiOutput(ns("ui_fst_out_status")),
            downloadButton(ns("run_FST_Analysis"), " Run",
                         icon  = icon("rocket"),
                         class = "btn-action-primary btn-block",
                         style = "font-weight: bold;"),
            run_folder_warning()
          ),
          column(9,
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fst_progress"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    )
  )
}