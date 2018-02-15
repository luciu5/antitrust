#library(shiny)
#library(rhandsontable)
require(shiny)
require(rhandsontable)

  ui <- shinyUI(fluidPage(

    titlePanel("Simulate a Merger"),
    sidebarLayout(
      sidebarPanel(
        helpText("Simulate a horizontal merger between the owners of products 1 and 2.",
                 "Enter shares, margins, and prices in Inputs table (right).",
                 "shares must be between 0 and 1.",
                 "Except for '2nd Score Auction', margins must also be be between 0 and 1."),

        checkboxInput("dispDetails", "Display Detailed Results", value = FALSE, width = NULL),
        checkboxInput("incEff", "Include Efficiencies", value = FALSE, width = NULL),
        #sliderInput("nProds", "# of Products:", min=3,max=20,value = 5,step=1),
        
        selectInput("supply", "Supply Specification:",
                    choices = c("Bertrand",
                                "2nd Score Auction",
                                "Cournot"
                                )),
        conditionalPanel(
          condition = "input.supply == 'Bertrand'",
          selectInput("demand_bert", "Demand Specification:",
                      choices = c("logit", "ces", "linear","aids"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Cournot'",
          selectInput("demand_cournot", "Demand Specification:",
                      choices = c("linear","log"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Cournot'",
          selectInput("cost", "Cost Specification:",
                      choices = c("linear","constant"))
        ),
        # conditionalPanel(
        #   condition = "input.supply == 'Cournot'",
        #   numericInput("mktPrice", "Market Price:",value=NA_real_, min=0)
        # ),
        conditionalPanel(
          condition = "input.supply == '2nd Score Auction'",
          selectInput("demand_2nd", "Demand Specification:",
                      choices = c("Logit"))
        )
      ),
      mainPanel(
        h2("Enter Inputs"),
        rHandsontableOutput("hot"),br(),
        actionButton("simulate","Simulate!"),
        br(), br(),
        h2("Results"),
        verbatimTextOutput("results")

      )
      
    )
  )
  )

  server <- shinyServer(function(input, output) {


   nPossProds <- 10 #only allow 10 products
   
   
   genInputData <- function(){

     
     inputData <- data.frame(
       Name = c("P1","P2","P3"),
       ownerPre  = c("O1","O2","O3"),
       ownerPost = c("O1","O1","O3"),
       Prices     = c(.0441,.0328,.0409),
       Shares    = c( 0.1344196, 0.3503055, 0.5152749),
       Margins   = c(.3830,.5515,.5421)
     )
     
     nDefProd <- nrow(inputData)
     inputData <- inputData[c(1:nDefProd,rep(1, nPossProds - nDefProd)),]
     if(input$incEff) inputData$mcDelta <- 0
     inputData[(nDefProd+1):nPossProds,] <- NA
     rownames(inputData) <- NULL
     
    
     return(inputData)

   }


   values <- reactiveValues(inputData = NULL)


    observe({

      
     values$inputData <- genInputData()

    })

    observe({

      input$supply

      if (!is.null(input$hot) && input$simulate == 0) {
        values$inputData = hot_to_r(input$hot)
      }
    })
    

    output$hot <- renderRHandsontable({
      inputData <- values[["inputData"]]
      if (!is.null(inputData))
        rhandsontable(inputData, stretchH = "all")
    })


    
    sims <- reactive({
      
      
      supply <- input$supply
      demand <- ifelse(supply =="Bertrand", input$demand_bert,
                       ifelse(supply == "Cournot", input$demand_cournot, input$demand_2nd))
      cost <- input$cost 
     
      indata <- values[["inputData"]]
      indata <- indata[!is.na(indata$Shares),]
      if(!input$incEff) indata$mcDelta <- 0
      
      indata$ownerPre <- factor(indata$ownerPre)
      indata$ownerPost <- factor(indata$ownerPost,levels=levels(indata$ownerPre))
      
      ownerPre = model.matrix(~-1+indata$ownerPre)
      ownerPre = tcrossprod(ownerPre)
      ownerPost = model.matrix(~-1+indata$ownerPost)
      ownerPost = tcrossprod(ownerPost)
      
    isolate({  
     
        switch(supply,
               Bertrand =
                 switch(demand,
                        logit= logit(prices= indata$Prices,
                                      shares= indata$Shares,
                                      margins= indata$Margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta),
                        ces= ces(prices= indata$Prices,
                                 shares= indata$Shares,
                                 margins= indata$Margins,
                                 ownerPre= ownerPre,
                                 ownerPost= ownerPost,
                                 mcDelta = indata$mcDelta),
                        linear=linear(prices= indata$Prices,
                                      quantities= indata$Shares,
                                      margins= indata$Margins,
                                      ownerPre= ownerPre,
                                      ownerPost= ownerPost,
                                      mcDelta = indata$mcDelta),
                        aids=aids(prices= indata$Prices,
                                  shares= indata$Shares,
                                  margins= indata$Margins,
                                  ownerPre= ownerPre,
                                  ownerPost= ownerPost,
                                  mcDelta = indata$mcDelta)
                 ),
               Cournot =                           cournot(prices= indata$Prices[1],
                                                           demand = demand,
                                                           cost= rep(cost, nrow(indata)),
                                                           quantities = as.matrix(indata$Shares),
                                                           margins= as.matrix(indata$Margins),
                                                           ownerPre= ownerPre,
                                                           ownerPost= ownerPost,
                                                           mcDelta = indata$mcDelta),
               `2nd Score Auction`=auction2nd.logit(prices= indata$Prices,
                                                    shares= indata$Shares,
                                                    margins= indata$Margins,
                                                    ownerPre= ownerPre,
                                                    ownerPost= ownerPost,
                                                    mcDelta = indata$mcDelta)
               
               
        )
      
      
    })
    })
    
    
    
    
    output$results <-
          
          renderPrint({
           
            #nQuant <- length(!is.na(values[["inputData"]]$Shares))
            #nQuant <- ifelse(length(nQuant)>0, nQuant,0)
            
            thisSim <- sims()
            if(input$simulate == 0) return(cat("Enter Inputs"))
            summary(thisSim,market= !input$dispDetails)
          })
        
      
    
    
    
  





  })

  ## run app
  shiny_antitrust <- function(){runApp(list(ui=ui, server=server))}


