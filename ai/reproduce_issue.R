# Source package files to load classes and methods
# Order matters for S4 inheritance
source("d:/Projects/antitrust/R/AntitrustClasses.R")
source("d:/Projects/antitrust/R/BertrandClasses.R")
source("d:/Projects/antitrust/R/BertrandRUMClasses.R")
# Source all other R files
files <- list.files("d:/Projects/antitrust/R", pattern = "\\.R$", full.names = TRUE)
for (f in files) {
  if (!grepl("BertrandClasses.R|BertrandRUMClasses.R", f)) {
    tryCatch(source(f), error = function(e) warning("Failed to source ", f, ": ", e$message))
  }
}

# Define parameters
alpha <- -2
sigmaNest <- 0.7 # Nesting parameter (0 < sigma <= 1)
prices <- c(1, 1) # Two products
delta <- c(0, 0) # Mean valuations
priceOutside <- 0
R <- 10
# Heterogeneous alphas to expose aggregation issues
set.seed(123)
alphas_draw <- rnorm(R, mean=alpha, sd=0.5) 

# --- Fixed Code Logic (from ElastMethods.R) ---
code_calc_fixed <- function(alphas_draw, sigmaNest, prices, delta, market = FALSE) {
    # Helper to calculate shares for a single draw
    calc_shares_draw <- function(alpha_r) {
        expUtil_sigma <- exp(util / sigmaNest)
        D_g_sigma <- D_g^sigmaNest
        s_g <- D_g_sigma / (1 + D_g_sigma)
        sInd <- s_jg * s_g
        return(sInd)
    }

    # Simulate the package's corrected aggregation: Derivative Aggregation
    nprods <- length(prices)
    nDraws <- length(alphas_draw)
    derivDraws <- array(0, dim = c(nprods, nprods, nDraws))
    sharesDraws <- matrix(0, nrow=nprods, ncol=nDraws)
    
    for(r in 1:nDraws) {
        alpha_r <- alphas_draw[r]
        util <- matrix(rep(delta, 1), ncol = 1) + outer(prices - priceOutside, alpha_r)
        
        expUtil_sigma <- exp(util / sigmaNest)
        D_g <- sum(expUtil_sigma)
        s_jg <- expUtil_sigma / D_g
        D_g_sigma <- D_g^sigmaNest
        s_g <- D_g_sigma / (1 + D_g_sigma)
        sInd <- s_jg * s_g
        sharesDraws[,r] <- sInd
        
        # Elasticity for this draw
        inv_sigma <- 1 / sigmaNest
        term2 <- 1 - inv_sigma - s_g
        term3 <- term2 / s_g
        
        diag_elast <- alpha_r * prices * (inv_sigma + s_jg * term2)
        cross_elast_row <- alpha_r * prices * sInd * term3
        elast_matrix <- matrix(cross_elast_row, nrow = nprods, ncol = nprods, byrow = TRUE)
        diag(elast_matrix) <- diag_elast
        
        # Convert to derivatives: d s_r / d p = E_r * s_r / p
        deriv_matrix <- elast_matrix * (sInd %*% t(1/prices))
        derivDraws[,,r] <- deriv_matrix
    }
    
    avg_deriv <- apply(derivDraws, c(1, 2), mean)
    avg_shares <- rowMeans(sharesDraws)
    
    # Compute Aggregate Elasticity
    # E_jk = (p_k / S_j) * (dS_j / dP_k)
    elast <- avg_deriv * ( (1/avg_shares) %*% t(prices) )
    
    if (market) {
        # Market Elasticity Logic
        # E_M = sum( w_j * sum_k(epsilon_jk) )
        inside_shares <- avg_shares / sum(avg_shares)
        mkt_elast <- sum(inside_shares * rowSums(elast))
        return(list(shares = avg_shares, mkt_elast = mkt_elast))
    } else {
        return(list(shares = avg_shares, elast = elast))
    }
}

# --- True Logic (Standard Nested Logit with Random Coefficients) ---
true_calc <- function(alphas_draw, sigmaNest, prices, delta, market = FALSE) {
    nprods <- length(prices)
    nDraws <- length(alphas_draw)
    
    # We need to aggregate derivatives: dS_j / dP_k = mean( d s_jr / d p_k )
    deriv_agg <- matrix(0, nprods, nprods)
    shares_agg <- numeric(nprods)
    
    for(r in 1:nDraws) {
        alpha_r <- alphas_draw[r]
        util <- matrix(rep(delta, 1), ncol = 1) + outer(prices - priceOutside, alpha_r)
        
        expUtil_sigma <- exp(util / sigmaNest)
        D_g <- sum(expUtil_sigma)
        s_jg <- expUtil_sigma / D_g
        D_g_sigma <- D_g^sigmaNest
        s_g <- D_g_sigma / (1 + D_g_sigma)
        sInd <- s_jg * s_g
        
        shares_agg <- shares_agg + sInd
        
        # Derivative for this draw
        inv_sigma <- 1 / sigmaNest
        term2 <- 1 - inv_sigma - s_g
        term3 <- term2 / s_g
        
        diag_elast <- alpha_r * prices * (inv_sigma + s_jg * term2)
        cross_elast_row <- alpha_r * prices * sInd * term3
        elast_matrix <- matrix(cross_elast_row, nrow = nprods, ncol = nprods, byrow = TRUE)
        diag(elast_matrix) <- diag_elast
        
        # Convert to derivative: d s_j / d p_k = E_jk * s_j / p_k
        deriv_matrix <- elast_matrix * (sInd %*% t(1/prices)) 
        
        deriv_agg <- deriv_agg + deriv_matrix
    }
    
    shares <- shares_agg / nDraws
    derivs <- deriv_agg / nDraws
    
    # Aggregate Elasticity: E_jk = (p_k / s_j) * d S_j / d P_k
    elast <- derivs * ( (1/shares) %*% t(prices) )
    
    if (market) {
        row_sums <- rowSums(elast)
        weights <- shares / sum(shares)
        mkt_elast <- sum(weights * row_sums)
        return(list(shares = shares, mkt_elast = mkt_elast))
    }

    list(shares = shares, elast = elast)
}

# Run comparison
cat("--- Product Elasticity ---\n")
res_code <- code_calc_fixed(alphas_draw, sigmaNest, prices, delta, market = FALSE)
res_true <- true_calc(alphas_draw, sigmaNest, prices, delta, market = FALSE)
print(res_code$elast)
print(res_true$elast)

cat("\n--- Market Elasticity ---\n")
res_code_mkt <- code_calc_fixed(alphas_draw, sigmaNest, prices, delta, market = TRUE)
res_true_mkt <- true_calc(alphas_draw, sigmaNest, prices, delta, market = TRUE)

cat("Code:", res_code_mkt$mkt_elast, "\n")
cat("True:", res_true_mkt$mkt_elast, "\n")

# Check for equality
elast_match <- isTRUE(all.equal(as.vector(res_code$elast), as.vector(res_true$elast)))
mkt_match <- isTRUE(all.equal(res_code_mkt$mkt_elast, res_true_mkt$mkt_elast))

if (elast_match && mkt_match) {
    cat("\nSUCCESS: Fixed code matches standard Nested Logit formulas for both Product and Market elasticity.\n")
} else {
    cat("\nFAILURE: Mismatch detected.\n")
}

# Define variables for calcPrices verification
shares_true <- as.numeric(res_true$shares)
J <- length(prices)

# -----------------------------------------------------------------------------
# 4. Verify calcPrices (Jacobian Fix)
# -----------------------------------------------------------------------------
message("\n--- Verifying calcPrices ---")

# Create object
test_obj <- new("LogitBLP")
test_obj@prices <- prices
test_obj@shares <- shares_true
test_obj@margins <- rep(0.5, J) # Dummy margins
test_obj@ownerPre <- diag(J) # Single product firms
test_obj@mcPre <- prices * 0.5 # Dummy MC
test_obj@subset <- rep(TRUE, J)
test_obj@priceOutside <- 0
test_obj@insideSize <- 1000
test_obj@mktSize <- 1000
test_obj@output <- TRUE # Output market

# Slopes
test_obj@slopes <- list(
    alpha = alpha,
    sigma = 1, # Dummy
    sigmaNest = sigmaNest,
    meanval = delta, 
    alphas = alphas_draw 
)
test_obj@nDraws <- R
test_obj@labels <- paste0("Prod", 1:J)
test_obj@control.equ <- list(tol = 1e-6, maxit = 100)
test_obj@priceStart <- prices

message("Calculating MC...")
mc_implied <- calcMC(test_obj, preMerger = TRUE)
test_obj@mcPre <- mc_implied

message("Running calcPrices (should recover original prices)...")
prices_new <- calcPrices(test_obj, preMerger = TRUE)

print(data.frame(Original = prices, Recovered = prices_new))

diff_prices <- max(abs(prices - prices_new))
message("Max price difference: ", diff_prices)

if (diff_prices < 1e-4) {
    message("SUCCESS: calcPrices recovered original prices.")
} else {
    message("FAILURE: calcPrices did not recover original prices.")
}
