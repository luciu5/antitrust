# Source package files
source("d:/Projects/antitrust/R/AntitrustClasses.R")
source("d:/Projects/antitrust/R/BertrandClasses.R")
source("d:/Projects/antitrust/R/BertrandRUMClasses.R")
files <- list.files("d:/Projects/antitrust/R", pattern = "\\.R$", full.names = TRUE)
for (f in files) {
  if (!grepl("BertrandClasses.R|BertrandRUMClasses.R|AntitrustClasses.R", f)) {
    tryCatch(source(f), error = function(e) warning("Failed to source ", f, ": ", e$message))
  }
}

# Parameters
alpha <- -2
sigmaNest <- 0.5 # Significant nesting
prices_pre <- c(1, 1)
prices_post <- c(1.1, 1.1) # Price increase
delta <- c(0, 0)
priceOutside <- 0
J <- 2
R <- 10
alphas_draw <- rep(alpha, R)

# Create LogitBLP object
test_obj <- new("LogitBLP")
test_obj@prices <- prices_pre
test_obj@pricePre <- prices_pre
test_obj@pricePost <- prices_post
test_obj@priceOutside <- priceOutside
test_obj@shares <- rep(0.1, J) # Dummy
test_obj@subset <- rep(TRUE, J)
test_obj@mktSize <- 1 # Normalized
test_obj@output <- TRUE
test_obj@slopes <- list(
  alpha = alpha,
  sigma = 1,
  sigmaNest = sigmaNest,
  meanval = delta,
  alphas = alphas_draw
)
test_obj@nDraws <- R

# Calculate CV using package
message("Calculating CV with package...")
cv_pkg <- CV(test_obj)
message("Package CV: ", cv_pkg)

# Calculate True CV
# CS = (1/alpha) * ln(1 + (sum exp(u/sigma))^sigma)
# VPre
utilPre <- matrix(rep(delta, R), nrow=R, byrow=TRUE) + outer(alphas_draw, prices_pre - priceOutside)
expUtilPre <- exp(utilPre / sigmaNest)
sumExpUtilPre <- rowSums(expUtilPre)
insideIVPre <- sumExpUtilPre^sigmaNest
VPre <- log(1 + insideIVPre)

# VPost
utilPost <- matrix(rep(delta, R), nrow=R, byrow=TRUE) + outer(alphas_draw, prices_post - priceOutside)
expUtilPost <- exp(utilPost / sigmaNest)
sumExpUtilPost <- rowSums(expUtilPost)
insideIVPost <- sumExpUtilPost^sigmaNest
VPost <- log(1 + insideIVPost)

# Individual CV
# CS_change = (VPost - VPre) / alpha
cv_ind_true <- (VPost - VPre) / alphas_draw
cv_true <- mean(cv_ind_true)

message("True CV: ", cv_true)

# Compare
diff_cv <- abs(cv_pkg - cv_true)
message("Difference: ", diff_cv)

if (diff_cv < 1e-6) {
  message("SUCCESS: Package CV matches standard formula.")
} else {
  message("FAILURE: Package CV does not match standard formula.")
  
  # Check if it matches with sigma multiplier
  cv_sigma <- mean(sigmaNest * (VPost - VPre) / alphas_draw)
  message("True CV * sigmaNest: ", cv_sigma)
  if (abs(cv_pkg - cv_sigma) < 1e-6) {
    message("CONFIRMED: Package CV has extra sigmaNest multiplier.")
  }
}
