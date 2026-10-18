# ui_import_data.R

import_data_ui <- function(id) {
  ns <- NS(id)

  box_title_style <- "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;"

  fluidPage(

    # CSS: equal-height top row boxes
    tags$style(HTML("
      .import-top-row > .row {
        display: flex !important;
        flex-wrap: wrap;
        align-items: stretch;
      }
      .import-top-row > .row > [class*='col-'] {
        display: flex !important;
        flex-direction: column;
      }
      .import-top-row > .row > [class*='col-'] > .box {
        flex: 1 1 auto !important;
      }
    ")),

    module_banner("upload", "Import Data",""),

    # ── Row 1: two config boxes side by side, same height ───────────────────
    tags$div(
      class = "import-top-row",
      fluidRow(

      column(
        width = 6,
        box(
          width = 12,
          title = div(style = box_title_style, icon("file-upload"), " Upload your data"),
          solidHeader = TRUE,

          fileInput(ns("file1"), "Choose CSV File", multiple = FALSE,
                    accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
          checkboxInput(ns("header"), "Header", TRUE),
          radioButtons(ns("sep"), "Separator",
                       choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                       selected = "\t"),
          tags$p(style = "color:#777;font-size:11px;", icon("info-circle"),
            " Choosing a file loads it automatically \u2014 no extra click needed."),
          uiOutput(ns("ui_load_status")),
          br(),
          actionButton(ns("load_default_data"), "Load Default Data", icon = icon("database"), class = "btn-action-secondary")
        )
      ),

      column(
        width = 6,
        box(
          width = 12,
          title = div(style = box_title_style, icon("sliders-h"), " Column assignment and formatting"),
          solidHeader = TRUE,
          footer = tagList("* mandatory fields", tags$br(), "Check column numbers below."),

          selectizeInput(ns("pop_data"),       "Population name*",    choices = NULL, options = list(placeholder = "select")),
          selectizeInput(ns("latitude_data"),  "Latitude",            choices = NULL, options = list(placeholder = "select")),
          selectizeInput(ns("longitude_data"), "Longitude",           choices = NULL, options = list(placeholder = "select")),
          tags$hr(),
          fluidRow(
            column(6,
              numericInput(ns("n_loci"), "Number of loci*", value = NA, min = 1, step = 1),
              tags$p(style = "color:#999;font-size:11px;margin-top:-6px;", "Auto-suggested once the first locus column is set below \u2014 adjust if needed.")),
            column(6, numericInput(ns("first_locus_col"), "Column number of the first locus*", value = NA, min = 1, step = 1))
          ),
          uiOutput(ns("ui_locus_range_preview")),
          textInput(ns("missing_code"), "Code for missing data", value = 0),
          br(),
          actionButton(ns("run_assign"), "Apply", icon = icon("check"), class = "btn-action-primary")
        )
      )
    ))  # closes fluidRow + tags$div(.import-top-row)
    ,
    fluidRow(
      box(
        width = 12,
        title = div(style = box_title_style, icon("list-ol"), " Column reference"),
        solidHeader = TRUE,
        tags$p(style = "color:#777;font-size:11px;",
          "Column numbers, for reference when filling in the fields above (e.g. the loci range)."),
        uiOutput(ns("ui_columns_index"))
      )
    )
  )
}