devtools::load_all()

prices <- c(10, 10)
shares <- c(0.3, 0.3)
margins <- c(0.8, 0.8)
ownerPre <- matrix(c(1,0,0,1), 2, 2)
ownerPost <- matrix(1, 2, 2)

cat("\n==========================================\n")
cat("         OUTPUT MARKET SIMULATIONS        \n")
cat("==========================================\n\n")

cat("--- LOGIT OUTPUT ---\n")
res_logit_out <- logit(prices, shares, margins, ownerPre=ownerPre, ownerPost=ownerPost, output=TRUE)
print(res_logit_out@slopes)
print(summary(res_logit_out))

cat("\n--- CES OUTPUT ---\n")
res_ces_out <- ces(prices, shares, margins, ownerPre=ownerPre, ownerPost=ownerPost, output=TRUE)
print(res_ces_out@slopes)
print(summary(res_ces_out))

cat("\n==========================================\n")
cat("         INPUT MARKET SIMULATIONS         \n")
cat("==========================================\n\n")

cat("--- LOGIT INPUT ---\n")
res_logit_in <- logit(prices, shares, margins, ownerPre=ownerPre, ownerPost=ownerPost, output=FALSE)
print(res_logit_in@slopes)
print(summary(res_logit_in))

cat("\n--- CES INPUT ---\n")
res_ces_in <- ces(prices, shares, margins, ownerPre=ownerPre, ownerPost=ownerPost, output=FALSE)
print(res_ces_in@slopes)
print(summary(res_ces_in))

cat("\n--- CES COURNOT INPUT ---\n")
res_ces_cournot_in <- ces.cournot(prices, shares, margins, ownerPre=ownerPre, ownerPost=ownerPost, output=FALSE)
print(res_ces_cournot_in@slopes)
print(summary(res_ces_cournot_in))

cat("\n--- AUCTION2ND CES INPUT ---\n")
res_auction_ces_in <- auction2nd.ces(prices, shares, margins, ownerPre=ownerPre, ownerPost=ownerPost, output=FALSE)
print(res_auction_ces_in@slopes)
print(summary(res_auction_ces_in))

cat("\n==========================================\n")
cat("      WELFARE COMPARISON (INPUT)          \n")
cat("==========================================\n\n")

cat("Logit Input   | CV (Supplier Harm): ", round(CV(res_logit_in), 4), "\n")
cat("CES   Input   | CV (Supplier Harm): ", round(CV(res_ces_in), 4), "\n")

cat("\nLogit Input   | PS (Buyer Profit): ", round(sum(calcProducerSurplus(res_logit_in, preMerger=FALSE)), 4), " (Post) vs ", round(sum(calcProducerSurplus(res_logit_in, preMerger=TRUE)), 4), " (Pre)\n")
cat("CES   Input   | PS (Buyer Profit): ", round(sum(calcProducerSurplus(res_ces_in, preMerger=FALSE)), 4), " (Post) vs ", round(sum(calcProducerSurplus(res_ces_in, preMerger=TRUE)), 4), " (Pre)\n")
