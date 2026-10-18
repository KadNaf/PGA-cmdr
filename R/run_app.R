#' Run the PGA-cmdr application
#'
#' @param ... Arguments passed to \code{\link[shiny]{shinyApp}}.
#' @return A \code{shiny.appobj} (invisibly).
#' @examples
#' \dontrun{
#'   run_app()
#' }
#' @export
run_app <- function(...) {
  options(shiny.maxRequestSize = 500 * 1024^2)  # 500 MB upload limit
  shiny::shinyApp(
    ui     = app_ui,
    server = app_server,
    ...
  )
}