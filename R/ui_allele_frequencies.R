# ui_allele_frequencies.R

ui_allele_frequencies <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    tags$head(gs_head()),
    
    module_banner("table", "Allele Frequencies",""),
    
    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("chart-pie"),
                    "Allele Frequency Analysis"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(4,
            h4(icon("save"), "Output file name"),
            tags$div(style = "max-width:320px;",
              textInput(ns("fstat_out_root"), NULL, value = "",
                        placeholder = "auto-filled from imported file")),
            uiOutput(ns("ui_fstat_out_status")),
            br(),
            downloadButton(ns("update_fstat"),
              label = "Run",
              icon = icon("rocket"),
              class = "btn-action-primary btn-block", style = "font-weight: bold;")
          )
        )
      )
    )
  )
}