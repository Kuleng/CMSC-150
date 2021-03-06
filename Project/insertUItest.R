
  # Define UI
  ui <- fluidPage(
    actionButton("add", "Add UI")
  )
  
  # Server logic
  server <- function(input, output, session) {
    observeEvent(input$add, {
      insertUI(
        selector = "#add",
        where = "afterEnd",
        ui = textInput(paste0("txt", input$add),
                       "Insert some text")
      )
    })
  }
  
  # Complete app with UI and server components
  shinyApp(ui, server)