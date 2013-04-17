rm(list=ls())
require(nleqslv)
require(ggplot2)
require(numDeriv)
currentdir <- getwd()
setwd("h:/TaragIC/AntiTrustRPackage/antitrust/R")
source("Antitrust.R")
source("Bertrand.R")
source("linear.R")
source("aids.R")
source("loglin.R")
source("logit.R")
source("logitALM.R")
source("logitNests.R")
source("logitNestsALM.R")
source("logitCap.R")
source("ces.R")
source("cesNests.R")
source("pcaids.R")
source("pcaidsNests.R")
source("sim.R")
source("auction2ndcap.R")
setwd(currentdir)


testMethods <- function(object,param){

    print(object)           # return predicted price change
    summary(object)         # summarize merger simulation

    print(elast(object,TRUE)  )    # returns premerger elasticities
    print(elast(object,FALSE) )    # returns postmerger elasticities

    print(diversion(object,TRUE))  # returns premerger diversion ratios
    print(diversion(object,FALSE)) # returns postmeger diversion ratios

    print(calcMargins(object,TRUE)) #margins
    print(calcMargins(object,FALSE))

    print(calcShares(object,TRUE)) #shares
    print(calcShares(object,FALSE))

    print(calcMC(object,TRUE)) #mc
    print(calcMC(object,FALSE))

    print(cmcr(object))            # returns the compensating marginal cost reduction


    print(plot(object))
    print(HypoMonTest(object,1:2)) # returns postmeger diversion ratios

    if(!missing(param)){
        print(CV(object,param))              # returns CV
    }

    else{ print(CV(object))  }            # returns CV





}


## Simulate a merger between two single-product firms in a
## three-firm market with loglin demand with diversions
## that are proportional to shares.
## This example assumes that the merger is between
## the first two firms



n <- 3 #number of firms in market
price    <- c(2.9,3.4,2.2)
quantity <- c(650,998,1801)
nests <- c("a","b","a")

slopes <- matrix(
c(-2.3,	  0.18,	0.28,
   0.11, -2.4,	0.1,
   0.13,   .16,-2.7),ncol=n)


margin <- -1/diag(slopes)


#simulate merger between firms 1 and 2
owner.pre <- diag(n)
owner.post <- owner.pre
owner.post[1,2] <- owner.post[2,1] <- 1

shares=quantity/sum(quantity)
## default is diversion according to revenue share
result1 <- aids(shares,margin,price,ownerPre=owner.pre,ownerPost=owner.post)
sim1 <- sim(price,demand="AIDS",list(slopes=result1@slopes,intercepts=result1@intercepts,mktElast=result1@mktElast),ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result1)
testMethods(sim1)

## User-supplied diversions (this will give the same answer as above)


result2 <- loglinear(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)
sim2 <- sim(price,demand="LogLin",list(slopes=result2@slopes,intercepts=result2@intercepts),ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result2)
testMethods(sim2)

## same as above, but with  pcaids demand
result3 <- pcaids(shares,-1/margin[1],-1,ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result3)

result4 <- linear(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)
sim4 <- sim(price,demand="Linear",list(slopes=result4@slopes,intercepts=result4@intercepts),ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result4)

##Beer calibration and simulation results from Epstein/Rubenfeld 2004, pg 80+
n <- 4 #number of firms in market
prodNames <- c("BUD","OLD STYLE","MILLER","MILLER-LITE","OTHER-LITE","OTHER-REG")
ownerPre <-c("BUD","OLD STYLE","MILLER","MILLER","OTHER-LITE","OTHER-REG")
ownerPost <-c("BUD","BUD","MILLER","MILLER","OTHER-LITE","OTHER-REG")
nests <- c("R","R","R","L","L","R")

price    <- c(.0441,.0328,.0409,.0396,.0387,.0497)
shares.quantity <- c(.066,.172,.253,.187,.099,.223)
shares.revenue <- c(.071,.137,.251,.179,.093,.269)
margins.logit <- c(.3830,.5515,.5421,.5557,.4453,.3769)
#margins.pcaids <- c(.4,.4179,.5208,.5208,.4059,.4589) #pcaids margins
margins.pcaids <- c(.4,.4232,.4997,.6787,.5725,.4787) #pcaids nests margins

names(price) <-
    names(shares.quantity) <-
    names(shares.revenue) <-
    names(margins.logit) <-
    names(margins.pcaids) <-
    prodNames


## ces demand
result5 <- ces(price,shares.revenue,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames,shareInside=.01)
demand.parm <- result5@slopes
demand.parm$shareInside=.01
sim5    <- sim(price,demand="CES",demand.parm,ownerPre=ownerPre,ownerPost=ownerPost)
testMethods(result5,1e9)

## Nested ces demand
result6 <- ces.nests(price,shares.revenue,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests,labels=prodNames,shareInside=.01)
result6.noconstraint <- ces.nests(price,shares.revenue,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests,labels=prodNames,shareInside=.01,constraint=FALSE)
demand.parm <- result6@slopes
demand.parm$shareInside=.01
sim6    <- sim(price,demand="CESNests",demand.parm,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests)
testMethods(result6,1e9)
testMethods(result6.noconstraint,1e9)

knownElast=-2.5
mktElast=-1.4


## PCAIDS demand
## Confirmed Against Epstein/Rubenfeld 2004 Beer example

result11 <- aids(shares.revenue,margins.pcaids,prices=price,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames)
testMethods(result11)

result7 <- pcaids(shares.revenue,knownElast,mktElast,knownElastIndex=1,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames)
testMethods(result7)

## Nested PCAIDS demand
## Confirmed Against Epstein/Rubenfeld 2004 Beer example
result8 <- pcaids.nests(shares.revenue,margins.pcaids,knownElast,mktElast,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests,labels=prodNames)
testMethods(result8)

## logit demand
result9 <- logit(price,shares.quantity,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames)
demand.parm <- result9@slopes
sim9    <- sim(price,demand="Logit",demand.parm,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames)
testMethods(result9)

## logit ALM demand
result12 <- logit.alm(price,shares.quantity,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames)
testMethods(result12)

## Nested logit demand
result10 <- logit.nests(price,shares.quantity,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests,labels=prodNames,parmsStart=c(-0.1117361,0.6340260))
result10.noconstraint <- logit.nests(price,shares.quantity,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests,labels=prodNames,constraint=FALSE)

demand.parm <- result10@slopes
sim10    <- sim(price,demand="LogitNests",demand.parm,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests,labels=prodNames)
testMethods(result10)

result10.2 <- logit.nests.alm(price,shares.quantity,margins.logit,ownerPre=ownerPre,ownerPost=ownerPost,nests=nests,labels=prodNames,parmsStart=c(-0.48110442,0.69168332,0.01972695))

## logit demand, with capacity constraints
mktSize=1000
cap <- c(66,200,300,200,99,300) # BUD and OTHER-LITE are capacity constrained
result11 <- logit.cap(price,shares.quantity,margins.logit,capacities=cap,mktSize=mktSize,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames)
cap <- c(200,200,300,200,99,300) # OTHER-LITE is capacity constrained
result11.2 <- logit.cap(price,shares.quantity,margins.logit,capacities=cap,mktSize=mktSize,ownerPre=ownerPre,ownerPost=ownerPost,labels=prodNames)
demand.parm <- result11@slopes
demand.parm$mktSize=mktSize
sim11    <- sim(price,demand="LogitCap",demand.parm,ownerPre=ownerPre,capacities=cap,ownerPost=ownerPost,labels=prodNames)


##
##Auctions
##



condShare=c(0.65,0.30,0.05)        #Vendor Shares, Conditional on Vendor Purchase
inShare=c(0.70)                    #Combined Vendor Share (equals 1 less In-House Share)
price=c(160,NA,NA)	           #Average Price of Each Vendor, in $000s
margin=c(0.65,NA,NA)           #Average Price-Cost Margin of Each Vendor


result.unif = auction2nd.cap(
          capacities=condShare,
          margins=margin,prices=price,reserve=NA,
          shareInside=inShare,
          sellerCostCDF="punif",constrain.reserve=FALSE,
          ownerPre=1:3,ownerPost=c(1,1,3)
          )

#Generate Fake data
result.unif@sellerCostParms <- result.unif@sellerCostBounds <- c(2,5)
result.unif@reservePre <- 4



shares=c(0.65,0.3,0.05)
shareIn.unif=cdfG(result.unif)
prices.unif=calcPrices(result.unif,exAnte=FALSE)
margin.unif=calcMargins(result.unif,exAnte=FALSE)



test.unif <- auction2nd.cap(
    capacities=shares,
    margins=margin.unif,prices=prices.unif,reserve=NA,
    shareInside=shareIn.unif,
    sellerCostCDF="punif",
    ownerPre=1:3,ownerPost=c(1,1,3)
    )

summary(test.unif)
str(test.unif)

test.exp <- auction2nd.cap(
    capacities=shares,
    margins=margin.unif,prices=prices.unif,
    shareInside=shareIn.unif,
    sellerCostCDF="pexp",reserve=NA,
    ownerPre=1:3,ownerPost=c(1,1,3)
    )


summary(test.exp)
str(test.exp)


## Weibull Distribution
test.w <- auction2nd.cap(
    capacities=condShare,
    margins=margin.unif,prices=prices.unif,
    shareInside=shareIn.unif,
    sellerCostCDF="pweibull",reserve=NA,
    ownerPre=1:3,ownerPost=c(1,1,3)
    )



summary(test.w)
str(test.w)


##Gumbel Distribution
test.g <- auction2nd.cap(
    capacities=shares,
    margins=margin.unif,prices=prices.unif,reserve=NA,
    shareInside=shareIn.unif,
    sellerCostCDF="pgumbel",constrain.reserve=FALSE,
    ownerPre=1:3,ownerPost=c(1,1,3)
    )



summary(test.g)
str(test.g)




##Frechet Distribution
test.f <- auction2nd.cap(
    capacities=shares,
    margins=margin.unif,prices=prices.unif,reserve=NA,
    shareInside=shareIn.unif,
    sellerCostCDF="pfrechet",parmsStart=c(4,1,1,3),
    ownerPre=1:3,ownerPost=c(1,1,3)
    )




summary(test.f)
str(test.f)
