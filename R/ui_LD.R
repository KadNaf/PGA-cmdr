# ui_LD.R

linkage_desequilibrium_UI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    useWaiter(),
    tags$head(
      tags$style(HTML("
        .ld-heatmap-container {
          overflow-x: auto;
          overflow-y: auto;
          max-height: 600px;
        }
      "))
    ),
    
    module_banner("link", "Linkage Disequilibrium",""),
    
    fluidRow(
      box(
        width = 4,
        title = div(
          style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
          icon("cogs"), "Linkage disequilibrium parameters"
        ),
        solidHeader = TRUE,
        status = "primary",
        
        h4(icon("sliders"), "Parameters"),
        numericInput(ns("n_iterations"), "Number of Permutations:",
                     value = 10000, min = 1000, max = 100000, step = 1000),
        
        tags$hr(),
        uiOutput(ns("ui_ld_out_status")),
        downloadButton(ns("run_LD"), " Run",
                     icon = icon("rocket"),
                     class = "btn-action-primary btn-block"),
            run_folder_warning()
      ),
      
      box(
        width = 8,
        title = div(
          style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
          icon("chart-line"), "Progress"
        ),
        solidHeader = TRUE,
        status = "primary",
        
        fluidRow(
          column(
            12,
            h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
            shinyWidgets::progressBar(id = ns("LD_progress"), value = 0, title = "Overall Progress")
          )
        ),
        
        br()
      )
    )
  )
}
