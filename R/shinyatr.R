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
          selectInput("demand", "Demand Specification:",
                      choices = c("logit", "ces", "linear","aids"))
        ),
        conditionalPanel(
          condition = "input.supply == 'Cournot'",
          selectInput("demand", "Demand Specification:",
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
          selectInput("demand", "Demand Specification:",
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
      
    ),
    hr(),
    print("Created by Charles Taragin (2018).The views expressed herein are entirely those of the authors and should not be purported to reflect those of the U.S. Department of Justice.  This merger simulator is released into the
          public domain without warranty of any kind, expressed or implied.")
  )
  )

  server <- shinyServer(function(input, output) {


   nPossProds <- 10 #only allow 10 products

   genInputData <- function(n
                            #,s
                            ){

     
     inputData <- data.frame(
       Name = c("BUD","OLD STYLE","MILLER"),
       ownerPre  = c("BUD","OLD STYLE","MILLER"),
       ownerPost = c("BUD","BUD","MILLER"),
       Prices     = c(.0441,.0328,.0409),
       Shares    = c( 0.1344196, 0.3503055, 0.5152749),
       Margins   = c(.3830,.5515,.5421)
     )
     
     nDefProd <- nrow(inputData)
     inputData <- inputData[c(1:nDefProd,rep(1, nPossProds - nDefProd)),]
     if(input$incEff) inputData$mcDelta <- 0
     inputData[(nDefProd+1):nPossProds,] <- NA
     rownames(inputData) <- NULL
     
     # inputData <- data.frame(Firm=1:input$nProds,
     #                         Prices=NA_real_,
     #                         Shares=NA_real_,
     #                         Margins=NA_real_,
     #                         ownerPre=factor(1:n),
     #                         ownerPost=factor(1:n),
     #                         mcDelta=0)

     #inputData$ownerPost[inputData$ownerPost==2]=1
     
     
     #if(s =="Cournot"){
     #   colnames(inputData)[colnames(inputData) == "Shares"] = "Quantities"
     # inputData$Prices <- NULL
     #}

     return(inputData)

   }


   values <- reactiveValues(inputData = NULL)


    observe({

      nProds <- input$nProds
      #supply <- input$supply
      
     values$inputData <- genInputData(nProds#, input$supply
                                      )

    })

    observe({



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
      demand <- input$demand
      cost <-   input$cost
      mktPrice <- input$mktPrice
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
                        logit = logit(prices= indata$Prices,
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
               Cournot =
                 switch(demand,
                        linear = 
                          cournot(prices = .5, quantities = matrix(c(100,10,50)), margins = matrix(c(.5,NA,NA)),ownerPre  =diag(3),ownerPost=matrix(c(1,1,0,1,1,0,0,0,1),ncol=3),mcDelta=rep(0,3))
                          
                          # cournot(prices= mktPrice,
                          #                demand = "linear",
                          #                cost= rep(cost, nrow(indata)),
                          #                quantities = as.matrix(indata$Shares),
                          #                margins= as.matrix(indata$Margins),
                          #                ownerPre= ownerPre,
                          #                ownerPost= ownerPost,
                          #                mcDelta = indata$mcDelta),
                        # log = cournot(prices= input$mktPrice,
                        #                        demand =  "log",
                        #                        cost= rep(cost, nrow(indata)),
                        #                      quantities = as.matrix(indata$shares),
                        #                      margins= as.matrix(indata$Margins),
                        #                      ownerPre= ownerPre,
                        #                      ownerPost= ownerPost,
                        #                      mcDelta = indata$mcDelta)
                        ),
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
            
            if(input$simulate == 0) return(cat("Enter Inputs"))
            summary(sims(),market= !inputdispDetails
                    )
          })
        
      
    
    
    
  





  })

  ## run app
  shiny_antitrust <- function(){runApp(list(ui=ui, server=server))}


