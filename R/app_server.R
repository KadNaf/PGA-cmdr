#' Application Server
#'
#' @param input,output,session Shiny server arguments.
#' @importFrom shiny reactiveValues
#' @export
app_server <- function(input, output, session) {
  rv <- shiny::reactiveValues(
    raw            = NULL,
    data           = NULL,
    formatted_data = NULL
  )

  server_import_data("import",        rv)
  server_allele_frequencies("allele", rv)
  server_general_stats("general_stats", rv)
  server_LD("ld",                     rv)
  server_null_alleles("null_alleles", rv)
  server_isolation_by_distance("ibd", rv)
}
