require(shiny)
require(rhandsontable)
require(antitrust)

shinyUI(fluidPage(

    titlePanel("Simulate a Merger"),
    sidebarLayout(
      sidebarPanel(
       h5(tags$b("Directions:")),
        helpText(tags$ul(
                 tags$li("Copy and paste (or enter) information into Inputs table (right) to simulate a merger between 'Firm1' and 'Firm2'"),
                 tags$li(helpText("See the",tags$a(href="https://CRAN.R-project.org/package=antitrust", "antitrust"),"R package vignette for more details about the models used here." ))
                 #tags$li("Shares must be between 0 and 1."),
                 #tags$li("Margins should exclude fixed costs.")
                 )
                 ),hr(),

        #checkboxInput("dispDetails", "Display Detailed Results", value = FALSE, width = NULL),
        #checkboxInput("calcElast", "Calibrate Market Elasticity", value = FALSE, width = NULL),
        radioButtons("calcElast", "Calibrate model parameters using:",
                     choices = c("market elasticity and 1 or more margins",
                                 "2 or more margins"
                     )),
        
        conditionalPanel(
          condition = "input.calcElast == 'market elasticity and 1 or more margins'",
          numericInput("enterElast", "Enter Market Elasticity:", value=-1,min=-Inf,max=0,step=.1#, width='75%'
                       )
        ),hr(),
        #checkboxInput("incEff", "Include Proportional Cost Changes (negative values imply cost reductions)", value = FALSE, width = NULL),
        
        radioButtons("supply", "Supply Specification:",
                    choices = c("Bertrand",
                                "2nd Score Auction",
                                "Cournot"
                                )),
        conditionalPanel(
          condition = "input.supply == 'Bertrand' && input.calcElast == '2 or more margins'",
          selectInput("demand_bert_alm", "Demand Specification:",
                      choices = c("logit (unknown elasticity)", "ces (unknown elasticity)", 
                                  #"linear",
                                  "aids (unknown elasticity)"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Bertrand'  && input.calcElast == 'market elasticity and 1 or more margins'",
          selectInput("demand_bert", "Demand Specification:",
                      choices = c("logit", "ces", 
                                  #"linear",
                                  "aids"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Cournot'  && input.calcElast == '2 or more margins' ",
          selectInput("demand_cournot_alm", "Demand Specification:",
                      choices = c("linear (unknown elasticity)","log (unknown elasticity)"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Cournot'  && input.calcElast  == 'market elasticity and 1 or more margins' ",
          selectInput("demand_cournot", "Demand Specification:",
                      choices = c("linear","log"))
        ),
        conditionalPanel(
          condition = "input.supply == '2nd Score Auction' && input.calcElast  ==  'market elasticity and 1 or more margins' ",
          selectInput("demand_2nd", "Demand Specification:",
                      choices = c("logit"))
        ),
        conditionalPanel(
          condition = "input.supply == '2nd Score Auction' && input.calcElast == '2 or more margins'",
          selectInput("demand_2nd_alm", "Demand Specification:",
                      choices = c("logit (unknown elasticity)"))
        ),hr(),
        conditionalPanel(
          condition = "input.supply == 'Cournot'",
          helpText(tags$b("Note:"), "only the first non-missing inputted price and product name is used for Cournot.")
          ),
        #conditionalPanel(
        #  condition = "input.supply != '2nd Score Auction'",
        #  helpText(tags$b("Note:"), "margins must be between 0 and 1.")
        #),
        #conditionalPanel(
        #  condition = "input.supply == '2nd Score Auction'",
        #  helpText(tags$b("Note:"), "2nd Score Auction requires $/unit margins.")
        #),
       conditionalPanel(
         condition = "input.supply == '2nd Score Auction' && input.calcElast  == 'market elasticity and 1 or more margins'",
         helpText(tags$b("Note:"), "2nd score Auction only requires a single price.")
       ),
       conditionalPanel(
         condition = "input.supply == '2nd Score Auction' && input.calcElast  == '2 or more margins'",
         helpText(tags$b("Note:"), "2nd score Auction does not require prices.")
       ),
        # conditionalPanel(
        #   condition = "input.demand_bert == 'aids'",
        #   helpText(tags$b("Note:"), "Only first non-missing inputted margin is used for pcaids. This margin should belong to a single product firm.")
        # ),
       conditionalPanel(
         condition = "input.demand_bert == 'aids'",
         helpText(tags$b("Note:"), "aids does not require pricing information.")
       )
      ),
      mainPanel(
        h2("Enter Inputs"),
        rHandsontableOutput("hot"), br(),
        #tags$head(
        #  tags$style(HTML('#run{color:white;background-color:black}'))
        #),
        actionButton(inputId ="simulate" , label = tags$b("Simulate"),style='padding:4px; font-size:125%')
      #)
        ,
        br(), br(),br(),
        tabsetPanel(id = "inTabset",
          tabPanel("Summary", value = "respanel", br(),br(),tableOutput("results"), br(),
                   helpText(tags$b("Note:"), "all price changes as well as compensating marginal cost reduction are (post-merger) share-weighted averages."),
                   conditionalPanel("input.demand_bert == 'aids' || input.demand_bert == 'ces' || input.demand_bert_alm == 'ces (unknown elasticity)'",
                                    helpText(tags$b("Note:"), "shares are revenue-based.")
                   )
          ),
          tabPanel("Details", value = "detpanel", br(),br(), tableOutput("results_shareOut"),br(), tableOutput("results_detailed"),
                   # conditionalPanel("input.demand_bert != 'aids' && input.supply != 'cournot'",
                   #    helpText(tags$b("Note:"), "shares include an outside good.")
                   #  ),
                   conditionalPanel("input.demand_bert == 'aids' || input.demand_bert == 'ces' || input.demand_bert_alm == 'ces (unknown elasticity)'",
                                    helpText(tags$b("Note:"), "shares are revenue-based.")
                   )), 
          tabPanel("Elasticities", value = "elastpanel",  br(),br(),
                   radioButtons("pre_elast", "",
                                                choices = c("Pre-Merger",
                                                            "Post-Merger"
                                                ), inline = TRUE),br(),
                   tableOutput("results_mktelast"),br(),
                   tableOutput("results_elast"),
                   conditionalPanel("input.supply != 'Cournot'",
                   helpText(tags$b("Note:"), "diagonal elements are own-price elasticities.","Off-diagonal elements are the cross-price elasticities of row with respect to column.")
                   ),
                   conditionalPanel("input.supply == 'Cournot'",
                                    helpText(tags$b("Note:"), "above are own-price elasticities")
                   )
                   ),
          tabPanel("Diagnostics", value = "diagpanel", br(),br(), h4("% Difference between predicted and observed values"), 
                   tableOutput("results_diag_elast"),
                   tableOutput("results_diagnostics"),
                   helpText(tags$b("Note:"), "Negative numbers mean that observed values are larger than predicted values."),br(),
                   h4("Parameters"),verbatimTextOutput("parameters"),
                   helpText("See the",tags$a(href="https://CRAN.R-project.org/package=antitrust", "antitrust"),"R package vignette for more details about the parameters displayed here." )
          ), 
          tabPanel("Messages", value = "msgpanel", br(),h4("Warnings"),  verbatimTextOutput("warnings"), br(),h4("Errors"),  verbatimTextOutput("errors"))
          
        )
        
      )
      
    )
  )
  )

 