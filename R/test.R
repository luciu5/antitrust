rm(list=ls())

currentdir <- getwd()
setwd("h:/TaragIC/AntiTrustRPackage/antitrust/R")
source("Classes.R")
source("Methods.R")
source("pclinear.R")
source("pcloglog.R")
source("linear.R")
source("loglog.R")
source("ces.R")
source("cesNests.R")
source("logit.R")
source("pcaids.R")
source("pcaidsNests.R")
setwd(currentdir)


testMethods <- function(object){

    print(object)           # return predicted price change
    summary(object)         # summarize merger simulation

    print(elast(object,TRUE)  )    # returns premerger elasticities
    print(elast(object,FALSE) )    # returns postmerger elasticities

    print(diversion(object,TRUE))  # returns premerger diversion ratios
    print(diversion(object,FALSE)) # returns postmeger diversion ratios

    print(cmcr(object))            # returns the compensating marginal cost reduction


}


## Simulate a merger between two single-product firms in a
## three-firm market with loglog demand with diversions
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
##result1 <- pcloglog(price,quantity,margin,shares=shares,ownerPre=owner.pre,ownerPost=owner.post)
##testMethods(result1)

## User-supplied diversions (this will give the same answer as above)



d=tcrossprod(1/(1-shares),shares)
diag(d)=1

result2 <- loglog(price,quantity,margin,diversions=d,ownerPre=owner.pre,ownerPost=owner.post)

testMethods(result2)


## same as above, but with linear demand
##result3 <- pclinear(price,quantity,margin,shares=shares,ownerPre=owner.pre,ownerPost=owner.post)
##testMethods(result3)

result4 <- linear(price,quantity,margin,diversions=d,ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result4)

## ces demand
result5 <- ces(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result5)

## Nested ces demand
result6 <- ces.nests(price,quantity,margin,ownerPre=owner.pre,ownerPost=owner.post,nests=nests)
testMethods(result6)


knownElast=-3
mktElast=-1

## PCAIDS demand
result8 <- pcaids(shares,knownElast,mktElast,ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result8)




price <- c(1,.9,.9,.9,0)
shares <- c(.808,.0822,.0402,.00374,.0658)

## logit demand
result7 <- logit(price,shares,margin,ownerPre=owner.pre,ownerPost=owner.post)
testMethods(result7)
