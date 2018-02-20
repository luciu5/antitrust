require(shiny)
require(rhandsontable)

shinyUI(fluidPage(

    titlePanel("Simulate a Merger"),
    sidebarLayout(
      sidebarPanel(
        helpText("Simulate a horizontal merger between the owners of products 1 and 2.",
                 "Enter (or copy and paste) shares, margins, and prices in Inputs table (right).",
                 "shares must be between 0 and 1.",
                 "Except for '2nd Score Auction', margins must also be be between 0 and 1."),

        checkboxInput("dispDetails", "Display Detailed Results", value = FALSE, width = NULL),
        #checkboxInput("incEff", "Include Proportional Cost Changes (negative values imply cost reductions)", value = FALSE, width = NULL),
        
        radioButtons("supply", "Supply Specification:",
                    choices = c("Bertrand",
                                "2nd Score Auction",
                                "Cournot"
                                )),
        conditionalPanel(
          condition = "input.supply == 'Bertrand'",
          selectInput("demand_bert", "Demand Specification:",
                      choices = c("logit", "ces", 
                                  #"linear",
                                  "aids"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Cournot'",
          selectInput("demand_cournot", "Demand Specification:",
                      choices = c("linear","log"))
        ),
       
        conditionalPanel(
          condition = "input.supply == '2nd Score Auction'",
          selectInput("demand_2nd", "Demand Specification:",
                      choices = c("Logit"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Cournot'",
          helpText("Note: only the first inputted price and product name is used for Cournot.")
          )
      ),
      mainPanel(
        h2("Enter Inputs"),
        rHandsontableOutput("hot"),
        #br(),
        #actionButton("simulate","Simulate!"),
        br(), br(),
        h2("Results"),
        tableOutput("results"),
        conditionalPanel(
          condition = "input.dispDetails == true",
          br(), br(),
          h2("Details"),
          tableOutput("results_detailed")
        )
      )
      
    )
  )
  )

 