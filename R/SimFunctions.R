#' @title Merger Simulation With User-Supplied Demand Parameters
#' @name Sim-Functions
#' @aliases sim
#' @description Simulates the price effects of a merger between two firms
#' with user-supplied demand parameters under the
#' assumption that all firms in the market are playing either a
#' differentiated products Bertrand pricing game, Cournot quantity game,
#' 2nd price (score) auction, or bargaining game.
#' @description Let k denote the number of products produced by all firms below.
#'
#' @param prices A length k vector of product prices.
#' @param shares A length k vector of product shares. Only used for \sQuote{BLP} demand. See Details.
#' @param supply A character string indicating how firms compete with one another. Valid
#' values are "bertrand" (Nash Bertrand), "cournot" (Nash Cournot), "auction2nd"
#' (2nd score auction), "bargaining", or "bargaining2nd".
#' @param demand A character string indicating the type of demand system
#'   to be used in the merger simulation. Supported demand systems are
#'   linear (\sQuote{Linear}), log-linear(\sQuote{LogLin}), logit (\sQuote{Logit}), nested logit
#'   (\sQuote{LogitNests}), ces (\sQuote{CES}), nested CES (\sQuote{CESNests}) and capacity
#'   constrained Logit (\sQuote{LogitCap}).
#' @param demand.param  See Below.
#' @param ownerPre EITHER a vector of length k whose values
#'   indicate which firm produced a product pre-merger OR
#'   a k x k matrix of pre-merger ownership shares.
#' @param ownerPost EITHER a vector of length k whose values
#'   indicate which firm produced a product after the merger OR
#'   a k x k matrix of post-merger ownership shares.
#' @param nests A length k vector identifying the nest that each
#'   product belongs to. Must be supplied when \sQuote{demand} equals \sQuote{CESNests} and
#'   \sQuote{LogitNests}.
#' @param capacities A length k vector of product capacities. Must be
#'   supplied when \sQuote{demand} equals \sQuote{LogitCap}.
#' @param mcDelta A vector of length k where each element equals the
#'   proportional change in a product's marginal costs due to
#'   the merger. Default is 0, which assumes that the merger does not
#'   affect any products' marginal cost.
#' @param subset A vector of length k where each element equals TRUE if
#'   the product indexed by that element should be included in the
#'   post-merger simulation and FALSE if it should be excluded.Default is a
#'   length k vector of TRUE.
#' @param insideSize A length 1 vector equal to total units sold if \sQuote{demand} equals "logit", or total revenues if
#' \sQuote{demand} equals "ces".
#' @param priceOutside A length 1 vector indicating the price of the
#'   outside good. This option only applies to the \sQuote{Logit} class and its child classes
#'   Default for \sQuote{Logit},\sQuote{LogitNests}, and \sQuote{LogitCap} is 0,
#'   and for \sQuote{CES} and \sQuote{CesNests} is 1.
#' @param priceStart A length k vector of starting values used to solve for
#'   equilibrium price. Default is the \sQuote{prices} vector for all values of
#'   demand except for \sQuote{AIDS}, which is set equal to a vector of 0s.
#' @param bargpowerPre A length k vector of pre-merger bargaining power parameters. Values
#' must be between 0 (sellers have the power) and 1 (buyers the power). Ignored if \sQuote{supply} not equal
#' to "bargaining" or bargaining2nd.
#' @param bargpowerPost A length k vector of post-merger bargaining power parameters. Values
#' must be between 0 (sellers have the power) and 1 (buyers the power). Default is \sQuote{bargpowerPre}.
#' Ignored if \sQuote{supply} not equal to "bargaining".
#' @param labels A k-length vector of labels. Default is \dQuote{Prod#}, where
#'   \sQuote{#} is a number between 1 and the length of \sQuote{prices}.
#' @param ... Additional options to feed to the
#'       optimizer used to solve for equilibrium prices.
#'
#' @details Using user-supplied demand parameters,
#' \code{sim} simulates the effects of a merger in a market where
#' firms are playing a differentiated products pricing game.
#'
#' The \sQuote{supply} parameter determines the type of competition.
#' When \sQuote{supply} equals \sQuote{cournot}, firms compete on quantities
#' rather than prices. Cournot competition is supported for \sQuote{Logit}
#' and \sQuote{BLP} demand systems, creating \code{LogitCournot} and
#' \code{CournotBLP} class objects respectively.
#'
#' If \sQuote{demand} equals \sQuote{Linear}, \sQuote{LogLin}, or
#' \sQuote{AIDS}, then \sQuote{demand.param} must be a
#' list containing:
#' \describe{
#'   \item{slopes}{A k x k matrix of slope coefficients.}
#'   \item{intercepts}{A length-k vector of intercepts.}
#' }
#' Additionally, if \sQuote{demand} equals \sQuote{AIDS},
#' \sQuote{demand.param} must also contain:
#' \describe{
#'   \item{mktElast}{An estimate of aggregate market elasticity.}
#' }
#' For \sQuote{Linear} demand models, \code{sim} returns an error if
#' any intercepts are negative. For \sQuote{Linear}, \sQuote{LogLin},
#' and \sQuote{AIDS} models, \code{sim} returns an error if not all
#' diagonal elements of the slopes matrix are negative.
#'
#' If \sQuote{demand} equals \sQuote{Logit} or \sQuote{LogitNests}, then
#' \sQuote{demand.param} must equal a list containing:
#' \describe{
#'   \item{alpha}{The price coefficient.}
#'   \item{meanval}{A length-k vector of mean valuations. If none of the
#'   values of \sQuote{meanval} are zero, an outside good is assumed to exist.}
#' }
#'
#' If \sQuote{demand} equals \sQuote{BLP}, then
#' \sQuote{demand.param} must equal a list containing:
#' \describe{
#'   \item{alpha}{The mean price coefficient (or use alphaMean).}
#'   \item{meanval}{A length-k vector of mean valuations. If none of the
#'   values of \sQuote{meanval} are zero, an outside good is assumed to exist.
#'   If \sQuote{meanval} is not provided, then \sQuote{shares} must be supplied
#'   so that mean valuations can be recovered via BLP contraction mapping.}
#'   \item{sigma}{The standard deviation of the random coefficient on price,
#'   representing consumer heterogeneity in price sensitivity.}
#'   \item{sigmaNest}{Optional nesting parameter for the outside good: sigmaNest in (0,1]
#'   where sigmaNest=1 is flat logit (no nesting) and sigmaNest->0 means products are
#'   perfect substitutes within the nest. Default is 1.}
#'   \item{piDemog}{Optional vector of demographic coefficients for the price coefficient.
#'   Each element represents the interaction effect of a demographic variable with price.}
#'   \item{demogMean}{Optional: Vector of length equal to piDemog length, containing the mean
#'   of each demographic variable. Default is 0 (demeaned demographics). If demographics were
#'   demeaned in estimation, use 0; otherwise provide the actual mean.}
#'   \item{demogCov}{Optional: Covariance matrix (square, dimensions = length of piDemog)
#'   for demographic variables. Default is identity matrix (unit variance, independent).
#'   Should match the variance structure of demographics in your data. For a single
#'   demographic with variance sigma^2, use matrix(sigma^2, nrow=1, ncol=1).}
#'   \item{nDraws}{Number of draws to use for simulating consumer heterogeneity. Default is 1000.}
#'   \item{prodChar}{Optional: k x L matrix of L product characteristics for k products.}
#'   \item{beta}{Optional: Length-L vector of mean coefficients on product characteristics.}
#'   \item{sigmaChar}{Optional: Length-L vector of random coefficient standard deviations on characteristics.}
#'   \item{pi}{Optional: nDemog x L matrix of demographic interactions with characteristics.}
#' }
#'
#' Note: The \sQuote{shares} argument is only used with \sQuote{BLP} demand.
#' If supplied for any other demand system, a warning will be issued but the function will proceed
#' (shares will be ignored). For \sQuote{BLP},
#' either \sQuote{meanval} must be provided in \sQuote{demand.param} OR \sQuote{shares}
#' must be supplied to the \code{sim} function to perform BLP contraction mapping.
#'
#' If \sQuote{demand} equals \sQuote{CES} or \sQuote{CESNests}, then
#' \sQuote{demand.param} must equal a list containing:
#' \describe{
#'   \item{gamma}{The price coefficient.}
#'   \item{alpha}{The coefficient on the numeraire good. May instead be
#'   calibrated using \sQuote{shareInside}.}
#'   \item{meanval}{A length-k vector of mean valuations. If none of the
#'   values of \sQuote{meanval} are zero, an outside good is assumed to exist.}
#'   \item{shareInside}{The budget share of all products in the
#'   market. Default is 1, meaning that all consumer wealth is spent on
#'   products in the market. May instead be specified using \sQuote{alpha}.}
#' }
#' @return \code{sim} returns an instance of the class specified by the
#' \sQuote{demand} argument.
#' @seealso The S4 class documentation for: \code{\linkS4class{Linear}},
#' \code{\linkS4class{AIDS}}, \code{\linkS4class{LogLin}}, \code{\linkS4class{Logit}},
#' \code{\linkS4class{LogitNests}}, \code{\linkS4class{CES}}, \code{\linkS4class{CESNests}}
#' @author Charles Taragin \email{ctaragin+antitrustr@gmail.com}
#'
#' @examples ## Calibration and simulation results from a merger between Budweiser and
#' ## Old Style. Note that the in the following model there is no outside
#' ## good; BUD's mean value has been normalized to zero.
#'
#' ## Source: Epstein/Rubenfeld 2004, pg 80
#'
#'
#' prodNames <- c("BUD", "OLD STYLE", "MILLER", "MILLER-LITE", "OTHER-LITE", "OTHER-REG")
#' ownerPre <- c("BUD", "OLD STYLE", "MILLER", "MILLER", "OTHER-LITE", "OTHER-REG")
#' ownerPost <- c("BUD", "BUD", "MILLER", "MILLER", "OTHER-LITE", "OTHER-REG")
#' nests <- c("Reg", "Reg", "Reg", "Light", "Light", "Reg")
#'
#' price <- c(.0441, .0328, .0409, .0396, .0387, .0497)
#'
#' demand.param <- list(
#'   alpha = -48.0457,
#'   meanval = c(0, 0.4149233, 1.1899885, 0.8252482, 0.1460183, 1.4865730)
#' )
#'
#' sim.logit <- sim(price,
#'   supply = "bertrand", demand = "Logit", demand.param = demand.param,
#'   ownerPre = ownerPre, ownerPost = ownerPost
#' )
#'
#'
#' print(sim.logit) # return predicted price change
#' summary(sim.logit) # summarize merger simulation
#'
#' elast(sim.logit, TRUE) # returns premerger elasticities
#' elast(sim.logit, FALSE) # returns postmerger elasticities
#'
#' diversion(sim.logit, TRUE) # return premerger diversion ratios
#' diversion(sim.logit, FALSE) # return postmerger diversion ratios
#'
#'
#' cmcr(sim.logit) # calculate compensating marginal cost reduction
#' upp(sim.logit) # calculate Upwards Pricing Pressure Index
#'
#' CV(sim.logit) # calculate representative agent compensating variation
#'
#' \dontrun{
#' ## Philadelphia National Bank Example: Comparing Four Model Types
#' ## Demonstrates Logit, Cournot, LogitBLP, and CournotBLP
#'
#' # Philadelphia National Bank data (38 products)
#' prices <- c(
#'   2.46, 2.21, 2.4, 2.29, 2.05, 2.45, 1.13, 2.39, 2.66, 2.73,
#'   1.83, 1.84, 2.76, 2.46, 1.1, 2.53, 2.27, 2.31, 2.72, 1.96, 2.26,
#'   2.63, 2.07, 2.31, 1.33, 2.67, 2.68, 2.59, 2.62, 1.85, 2.63, 1.51,
#'   1.8, 2.66, 1.4, 2.39, 2.35, 1.82
#' )
#'
#' shares <- c(
#'   0.12499, 0.11462, 0.07748, 0.05231, 0.05168, 0.02918, 0.01622,
#'   0.01245, 0.00992, 0.00578, 0.00461, 0.00417, 0.00415, 0.00275,
#'   0.00271, 0.00214, 0.00206, 0.00184, 0.00176, 0.00155, 0.00151,
#'   0.00146, 0.00146, 0.00103, 0.00101, 0.00095, 8e-04, 0.00066,
#'   0.00062, 0.00061, 0.00057, 0.00053, 0.00053, 5e-04, 0.00049,
#'   0.00048, 0.00041, 0.00039
#' )
#'
#' ownerPre <- diag(38)
#' ownerPost <- diag(38)
#' ownerPost[2, 3] <- 1 # Products 2 and 3 merge
#' ownerPost[3, 2] <- 1
#'
#' insideSize <- 4773473000
#'
#' # Logit parameters (from calibrated model)
#' logit_params <- list(
#'   alpha = 0.4296,
#'   meanval = c(
#'     -2.265, -2.206, -2.671, -2.972, -2.855, -3.596, -3.494, -4.344,
#'     -4.709, -5.253, -4.966, -5.059, -5.572, -5.773, -5.128, -6.043,
#'     -5.933, -6.057, -6.339, -6.027, -6.215, -6.46, -6.146, -6.597,
#'     -6.131, -6.871, -7.035, -7.169, -7.237, -6.833, -7.318, -6.804,
#'     -6.939, -7.457, -6.819, -7.328, -7.473, -7.227
#'   )
#' )
#'
#' # BLP parameters
#' blp_params <- list(
#'   alpha = 1.0,
#'   sigma = 0.3,
#'   sigmaNest = 0.9,
#'   nDraws = 200
#' )
#'
#' # MODEL 1: Logit (Bertrand)
#' result_logit <- sim(prices,
#'   supply = "bertrand", demand = "Logit",
#'   demand.param = logit_params, ownerPre = ownerPre,
#'   ownerPost = ownerPost, insideSize = insideSize
#' )
#' summary(result_logit)
#'
#' # MODEL 2: Logit (Cournot)
#' result_cournot <- sim(prices,
#'   supply = "cournot", demand = "Logit",
#'   demand.param = logit_params, ownerPre = ownerPre,
#'   ownerPost = ownerPost, insideSize = insideSize
#' )
#' summary(result_cournot)
#'
#' # MODEL 3: LogitBLP (Bertrand with random coefficients)
#' result_logitblp <- sim(prices, shares,
#'   supply = "bertrand", demand = "BLP",
#'   demand.param = blp_params, ownerPre = ownerPre,
#'   ownerPost = ownerPost, insideSize = insideSize
#' )
#' summary(result_logitblp)
#'
#' # MODEL 4: CournotBLP (Cournot with random coefficients)
#' result_cournotblp <- sim(prices, shares,
#'   supply = "cournot", demand = "BLP",
#'   demand.param = blp_params, ownerPre = ownerPre,
#'   ownerPost = ownerPost, insideSize = insideSize
#' )
#' summary(result_cournotblp)
#' }
#'
#' @include LogitFunctions.R
NULL


#' @rdname Sim-Functions
#' @export
sim <- function(prices,
                shares = NULL,
                supply = c("bertrand", "cournot", "auction", "bargaining", "bargaining2nd"),
                demand = c("Linear", "AIDS", "LogLin", "Logit", "CES", "LogitNests", "CESNests", "LogitCap", "BLP"), demand.param,
                ownerPre, ownerPost, nests, capacities,
                mcDelta = rep(0, length(prices)),
                subset = rep(TRUE, length(prices)),
                insideSize = 1,
                priceOutside,
                priceStart,
                bargpowerPre = rep(0.5, length(prices)),
                bargpowerPost = bargpowerPre,
                labels = paste("Prod", 1:length(prices), sep = ""),
                ...) {
  demand <- match.arg(demand)
  supply <- match.arg(supply)
  nprods <- length(prices)

  # Validate supply/demand combinations
  valid_combinations <- list(
    bertrand = c("Linear", "AIDS", "LogLin", "Logit", "CES", "LogitNests", "CESNests", "LogitCap", "BLP"),
    cournot = c("Logit", "CES", "BLP"),
    auction = c("Logit", "CES"),
    bargaining = c("Logit"),
    bargaining2nd = c("Logit")
  )

  if (!(demand %in% valid_combinations[[supply]])) {
    stop(paste(supply, "/", demand, "currently not supported."))
  }


  if (missing(priceStart)) {
    if (demand == "AIDS") {
      priceStart <- runif(nprods)
    } else {
      priceStart <- prices
    }
  }

  ## Create placeholders values to fill required Class slots

  sharesProvided <- !is.null(shares) && !missing(shares)

  if (is.null(shares) || missing(shares)) {
    shares <- rep(1 / nprods, nprods)
  } else {
    if (length(shares) != nprods) {
      stop("'shares' must have the same length as 'prices'")
    }
    # Warn if shares provided for non-BLP demand
    if (demand != "BLP") {
      warning("'shares' argument is only used for BLP demand (for BLP contraction mapping). For '", demand, "', shares will be ignored.")
    }
  }

  margins <- rep(1 / nprods, nprods)


  if (!missing(nests)) {
    nests <- factor(nests, levels = unique(nests))
  }


  ## general checks
  if (!is.list(demand.param)) {
    stop("'demand.param' must be a list.")
  }

  ## Checks for discrete choice models
  if (demand %in% c("CESNests", "LogitNests", "CES", "Logit", "LogitCap", "BLP")) {
    # meanval is optional for BLP (recovered via contraction if not provided)
    if (demand != "BLP") {
      if (!("meanval" %in% names(demand.param))) {
        stop("'demand.param' does not contain 'meanval'.")
      }
      if (length(demand.param$meanval) != nprods || any(is.na(demand.param$meanval))) {
        stop("'meanval' must be a length-k vector of product mean valuations. NAs not allowed.")
      }
    } else {
      # For BLP, validate meanval only if provided
      if ("meanval" %in% names(demand.param)) {
        if (length(demand.param$meanval) != nprods || any(is.na(demand.param$meanval))) {
          stop("'meanval' must be a length-k vector of product mean valuations. NAs not allowed.")
        }
      }
    }

    if (demand %in% c("LogitNests", "Logit", "LogitCap", "BLP")) {
      if (demand == "BLP") {
        # Check if meanval is provided (optional for BLP)
        # If NOT provided, calcSlopes will recover it via BLP contraction from observed shares
        if (!("meanval" %in% names(demand.param))) {
          # If meanval not provided, shares MUST be provided for BLP contraction
          if (!sharesProvided) {
            stop("For BLP, either 'meanval' must be in 'demand.param' OR 'shares' must be provided. Cannot recover delta without observed shares.")
          }
          message("Note: 'meanval' (delta) not provided for BLP. It will be recovered via BLP contraction from observed shares/prices.")
        }

        # No need to check random_draws as they are handled in calcSlopes
        if (!("sigma" %in% names(demand.param)) || length(demand.param$sigma) != 1) {
          stop("'demand.param' must contain 'sigma', a scalar parameter for price coefficient heterogeneity.")
        }
        # Set default nDraws if not provided
        if (!("nDraws" %in% names(demand.param))) {
          demand.param$nDraws <- 500
          message("'nDraws' not provided for BLP. Defaulting to 500 draws.")
        }
        # Handle demographic parameters
        # Infer nDemog from piDemog length (nDemog is redundant)
        if (!("piDemog" %in% names(demand.param))) {
          demand.param$piDemog <- numeric(0)
        }

        # Auto-detect nDemog from piDemog length (always inferred, never user-specified)
        demand.param$nDemog <- length(demand.param$piDemog)

        # Set defaults for demogMean and demogCov if not provided
        if (demand.param$nDemog > 0) {
          if (!("demogMean" %in% names(demand.param))) {
            demand.param$demogMean <- rep(0, demand.param$nDemog) # Default: demeaned
            message("demogMean not provided. Defaulting to 0 (demeaned demographic).")
          }
          if (!("demogCov" %in% names(demand.param))) {
            demand.param$demogCov <- diag(demand.param$nDemog) # Default: unit variance, independent
            warning(
              "demogCov not provided. Defaulting to identity matrix (unit variance). ",
              "Consider specifying demogCov based on your data's variance."
            )
          }
        }

        # Validate demographic parameters
        if ("demogMean" %in% names(demand.param)) {
          if (length(demand.param$demogMean) != demand.param$nDemog) {
            stop(
              "demogMean length (", length(demand.param$demogMean),
              ") does not match piDemog length (", demand.param$nDemog, ")."
            )
          }
        }
        if ("demogCov" %in% names(demand.param)) {
          if (!is.matrix(demand.param$demogCov) ||
            nrow(demand.param$demogCov) != demand.param$nDemog ||
            ncol(demand.param$demogCov) != demand.param$nDemog) {
            # Validate dimensions match
            if (nrow(demand.param$demogCov) != ncol(demand.param$demogCov)) {
              stop("demogCov must be a square matrix.")
            }
            # Check that demogMean and demogCov dimensions are consistent
            if (length(demand.param$demogMean) != nrow(demand.param$demogCov)) {
              stop(
                "demogMean length (", length(demand.param$demogMean),
                ") must match demogCov dimensions (", nrow(demand.param$demogCov), "x",
                ncol(demand.param$demogCov), ")."
              )
            }
            stop("demogCov must be a square matrix with dimensions equal to piDemog length.")
          }
        }
        if (!is.numeric(demand.param$sigmaNest) || length(demand.param$sigmaNest) != 1 ||
          is.na(demand.param$sigmaNest) || !is.finite(demand.param$sigmaNest) ||
          demand.param$sigmaNest <= 0 || demand.param$sigmaNest > 1) {
          stop("'sigmaNest' (nesting parameter) must be a single numeric value in (0,1].")
        }

        # NEW: Validate product characteristics parameters
        if ("prodChar" %in% names(demand.param)) {
          if (!is.matrix(demand.param$prodChar) || nrow(demand.param$prodChar) != nprods) {
            stop("'prodChar' must be a k x L matrix where k is the number of products.")
          }
          nChar <- ncol(demand.param$prodChar)

          if ("beta" %in% names(demand.param)) {
            if (length(demand.param$beta) != nChar) {
              stop("'beta' must have length equal to number of characteristics in 'prodChar'.")
            }
          } else {
            stop("If 'prodChar' is provided, 'beta' (mean coefficients) must also be provided.")
          }

          if ("sigmaChar" %in% names(demand.param)) {
            if (length(demand.param$sigmaChar) != nChar) {
              stop("'sigmaChar' must have length equal to number of characteristics.")
            }
          }

          if ("pi" %in% names(demand.param)) {
            if (demand.param$nDemog == 0) {
              warning("'pi' provided but 'nDemog' is 0. Ignoring 'pi'.")
              demand.param$pi <- NULL
            } else {
              if (!is.matrix(demand.param$pi) ||
                nrow(demand.param$pi) != demand.param$nDemog ||
                ncol(demand.param$pi) != nChar) {
                stop("'pi' must be a nDemog x L matrix where L is the number of characteristics.")
              }
            }
          }
        }
      }

      ## An outside option is assumed to exist if all mean valuations are non-zero
      ## For BLP, meanval might not be in demand.param yet (recovered via BLP contraction)
      if ("meanval" %in% names(demand.param)) {
        if (all(demand.param$meanval != 0)) {
          normIndex <- NA
          # Only override shares if they weren't explicitly provided by user
          if (!sharesProvided) {
            shares <- rep(1 / (nprods + 1), nprods)
          }
        } else {
          normIndex <- which(demand.param$meanval == 0)

          if (length(normIndex) > 1) {
            warning("multiple values of meanval are equal to zero. Normalizing with respect to the first product with zero mean value.")
            normIndex <- normIndex[1]
          }
        }
      } else {
        # meanval not provided - this should only happen for BLP
        # and we've already verified shares were provided (error thrown at line 251 if not)
        # Assume outside good exists, will be determined after BLP contraction
        normIndex <- NA
        # sharesProvided must be TRUE here (or we would have errored), so don't override shares
      }

      # Determine alpha mean and set output sign accordingly (TRUE if alpha<0, FALSE if alpha>0)
      # Accept aliases: alpha, alphaMean, alpha_mean
      alphaVal <- ifelse("alpha" %in% names(demand.param), demand.param$alpha,
        ifelse("alphaMean" %in% names(demand.param), demand.param$alphaMean,
          ifelse("alpha_mean" %in% names(demand.param), demand.param$alpha_mean, NA)
        )
      )

      if (length(alphaVal) != 1 || is.na(alphaVal)) {
        stop("'demand.param' must include a scalar 'alpha' (or 'alphaMean'/'alpha_mean').")
      }
      # Normalize into 'alpha' for downstream methods
      demand.param$alpha <- alphaVal
      # output=TRUE for negative alpha (output market), FALSE for positive alpha (input market)
      outputFlag <- alphaVal < 0

      shareInside <- sum(shares)
      if (missing(priceOutside)) {
        priceOutside <- 0
      }
    } else if (demand %in% c("CESNests", "CES")) {
      if (!("gamma" %in% names(demand.param)) ||
        length(demand.param$gamma) != 1 ||
        isTRUE(demand.param$gamma < 0)) {
        stop("'demand.param' does not contain 'gamma' or 'gamma' is not a positive number.")
      }


      ## uncover Numeraire Coefficients
      if (!("alpha" %in% names(demand.param)) &&
        !("shareInside" %in% names(demand.param))) {
        warning("'demand.param' does not contain either 'alpha' or 'shareInside'. Setting shareInside=1 and alpha=NULL.")
        shareInside <- 1
        demand.param$alpha <- NULL
      } else if ("shareInside" %in% names(demand.param)) {
        shareInside <- demand.param$shareInside
        demand.param$shareInside <- NULL

        if (shareInside < 1) {
          demand.param$alpha <- 1 / shareInside - 1
        } else {
          demand.param$alpha <- NULL
        }
      } else {
        shareInside <- 1 / (1 + demand.param$alpha)
      }

      outputFlag <- TRUE


      ## An outside option is assumed to exist if all mean valuations are non-zero
      if (all(demand.param$meanval != 1)) {
        normIndex <- NA
        shares <- rep(1 / (nprods + 1), nprods)
      } else {
        normIndex <- which(demand.param$meanval == 1)

        if (length(normIndex) > 1) {
          warning("multiple values of meanval are equal to one. Normalizing with respect to the first product with  mean value equal to 1.")
          normIndex <- normIndex[1]
        }
      }


      if (missing(priceOutside)) {
        priceOutside <- 1
      }
    }

    if (demand %in% c("CESNests", "LogitNests")) {
      if (!("sigma" %in% names(demand.param))) {
        stop("'demand.param' does not contain 'sigma'.")
      }

      if (missing(nests) ||
        length(nests) != nprods) {
        stop("When 'demand' equals 'CESNests' or 'LogitNests', 'nests' must equal a vector whose length equals the number of products.")
      }

      if (length(demand.param$sigma) == 1) {
        constraint <- TRUE
        demand.param$sigma <- rep(demand.param$sigma, nlevels(nests))
      } else {
        constraint <- FALSE
      }


      if (nlevels(nests) != length(demand.param$sigma)) {
        stop("The number of nests in 'nests' must either equal the number of nesting parameters in 'demand.param$sigma'.")
      }
    }
  }


  ## Checks for Linear-demand style models
  if (demand %in% c("Linear", "LogLin", "AIDS")) {
    if (!("slopes" %in% names(demand.param))) {
      stop("'demand.param' does not contain 'slopes'")
    }
    if (!("intercepts" %in% names(demand.param))) {
      stop("'demand.param' does not contain 'intercepts'")
    }

    if (!(is.matrix(demand.param$slopes)) ||
      ncol(demand.param$slopes) != nprods ||
      nrow(demand.param$slopes) != nprods ||
      any(diag(demand.param$slopes) > 0)) {
      stop("'slopes' must be a k x k matrix of slope coeficients whose diagonal elements must all be negative.")
    }
    if (!is.vector(demand.param$intercepts) ||
      length(demand.param$intercepts) != nprods ||
      isTRUE(any(demand.param$intercepts < 0, na.rm = TRUE))) {
      stop("'intercepts' must be a length-k vector whose elements are all non-negative")
    }

    if (demand == "AIDS" &&
      !("mktElast" %in% names(demand.param))) {
      warning("'demand.param' does not contain 'mktElast'. Setting 'mktElast' equal to -1")
      demand.param$mktElast <- -1
    }
  }


  ## Create constructors for each demand system specified in the 'demand' parameter

  if (demand == "CESNests") {
    result <- new(demand,
      prices = prices, shares = shares, margins = margins,
      mcDelta = mcDelta,
      subset = subset,
      ownerPre = ownerPre,
      ownerPost = ownerPost,
      nests = nests,
      normIndex = normIndex,
      parmsStart = c(demand.param$gamma, demand.param$sigma),
      priceStart = priceStart,
      constraint = constraint,
      shareInside = shareInside, labels = labels
    )
  } else if (demand == "LogitNests") {
    result <- new(demand,
      prices = prices, shares = shares, margins = margins,
      mcDelta = mcDelta,
      subset = subset,
      ownerPre = ownerPre,
      ownerPost = ownerPost,
      nests = nests,
      normIndex = normIndex,
      parmsStart = c(demand.param$alpha, demand.param$sigma),
      priceStart = priceStart,
      constraint = constraint,
      shareInside = shareInside, labels = labels
    )
  } else if (demand == "BLP") {
    # Determine class based on supply type
    if (supply == "cournot") {
      demandClass <- "CournotBLP"
    } else {
      demandClass <- "LogitBLP"
    }

    result <- new(demandClass,
      prices = prices,
      shares = shares,
      margins = margins,
      slopes = demand.param,
      normIndex = normIndex,
      mcDelta = mcDelta,
      insideSize = insideSize,
      subset = subset,
      ownerPre = ownerPre,
      ownerPost = ownerPost,
      priceStart = priceStart,
      priceOutside = priceOutside,
      shareInside = shareInside,
      output = outputFlag,
      nDraws = demand.param$nDraws,
      labels = labels
    )
  } else if (demand %in% c("Logit", "CES")) {
    result <- switch(supply,
      bertrand = new(demand,
        prices = prices, shares = shares,
        margins = margins,
        normIndex = normIndex,
        mcDelta = mcDelta,
        insideSize = insideSize,
        subset = subset,
        ownerPre = ownerPre,
        ownerPost = ownerPost,
        priceStart = priceStart,
        priceOutside = priceOutside,
        shareInside = shareInside,
        output = outputFlag,
        labels = labels
      ),
      cournot = new(paste0(demand, "Cournot"),
        prices = prices, shares = shares,
        margins = margins,
        normIndex = normIndex,
        mcDelta = mcDelta,
        insideSize = insideSize,
        subset = subset,
        ownerPre = ownerPre,
        ownerPost = ownerPost,
        priceStart = priceStart,
        priceOutside = priceOutside,
        shareInside = shareInside,
        output = outputFlag,
        labels = labels
      ),
      bargaining = new("BargainingLogit",
        prices = prices, shares = shares,
        margins = margins,
        normIndex = normIndex,
        ownerPre = ownerPre,
        ownerPost = ownerPost,
        bargpowerPre = bargpowerPre,
        bargpowerPost = bargpowerPost,
        insideSize = insideSize,
        mcDelta = mcDelta,
        subset = subset,
        priceOutside = priceOutside,
        shareInside = shareInside,
        priceStart = priceStart,
        output = outputFlag,
        labels = labels,
        cls = "BargainingLogit"
      ),
      auction = new(paste0("Auction2nd", demand),
        prices = prices, shares = shares,
        margins = margins,
        normIndex = normIndex,
        ownerPre = ownerPre,
        ownerPost = ownerPost,
        insideSize = insideSize,
        mcDelta = mcDelta,
        subset = subset,
        priceOutside = priceOutside,
        shareInside = shareInside,
        priceStart = priceStart,
        output = outputFlag,
        labels = labels,
        cls = paste0("Auction2nd", demand)
      ),
      bargaining2nd = new("Bargaining2ndLogit",
        prices = prices, shares = shares,
        margins = margins,
        normIndex = normIndex,
        ownerPre = ownerPre,
        ownerPost = ownerPost,
        bargpowerPre = bargpowerPre,
        bargpowerPost = bargpowerPost,
        insideSize = insideSize,
        mcDelta = mcDelta,
        subset = subset,
        priceOutside = priceOutside,
        shareInside = shareInside,
        priceStart = priceStart,
        output = outputFlag,
        labels = labels,
        cls = "Bargaining2ndLogit"
      )
    )
  } else if (demand == "LogitCap") {
    if (!("mktSize" %in% names(demand.param))) {
      if (!missing(capacities)) {
        warning("'demand.param' does not contain 'mktSize'. Setting 'mktSize' equal to the sum of 'capacities'.")
        mktSize <- sum(capacities)
      } else {
        stop("'demand.param' does not contain 'mktSize'")
      }
    } else {
      mktSize <- demand.param$mktSize
    }


    shares <- capacities / mktSize
    shares <- shares / sum(shares)

    result <- new(demand,
      prices = prices, shares = shares,
      margins = margins, capacities = capacities, mktSize = mktSize,
      normIndex = normIndex,
      ownerPre = ownerPre,
      ownerPost = ownerPost,
      mcDelta = mcDelta,
      subset = subset,
      priceStart = priceStart, shareInside = shareInside,
      output = outputFlag,
      labels = labels
    )
  } else if (demand == "Linear") {
    result <- new(demand,
      prices = prices, quantities = shares, margins = margins,
      shares = shares, mcDelta = mcDelta, subset = subset,
      ownerPre = ownerPre, diversion = -diag(nprods),
      symmetry = identical(demand.param$slopes, t(demand.param$slopes)),
      ownerPost = ownerPost, priceStart = priceStart, labels = labels
    )
  } else if (demand == "AIDS") {
    ## find the market elasticity that best explains user-supplied intercepts and prices

    aidsShares <- as.vector(demand.param$intercepts + demand.param$slopes %*% log(prices)) # AIDS needs actual shares for prediction
    aidsDiv <- tcrossprod(1 / (1 - aidsShares), aidsShares)
    diag(aidsDiv) <- -1

    result <- new(demand,
      prices = prices, quantities = shares, margins = margins,
      shares = aidsShares,
      mcDelta = mcDelta, subset = subset, mktElast = demand.param$mktElast,
      ownerPre = ownerPre, diversion = aidsDiv,
      priceStart = priceStart,
      ownerPost = ownerPost, labels = labels
    )
  } else if (demand == "LogLin") {
    result <- new(demand,
      prices = prices, quantities = shares, margins = margins,
      shares = shares, mcDelta = mcDelta, subset = subset, priceStart = priceStart,
      ownerPre = ownerPre, diversion = -diag(nprods),
      ownerPost = ownerPost, labels = labels
    )
  }


  if (demand %in% c("Linear", "LogLin", "AIDS")) {
    result@slopes <- demand.param$slopes
    result@intercepts <- demand.param$intercepts
  } else {
    result@slopes <- demand.param
  }

  ## For BLP, recover delta and generate random coefficients via calcSlopes
  if (demand == "BLP") {
    result <- calcSlopes(result)
  }


  ## Convert ownership vectors to ownership matrices
  result@ownerPre <- ownerToMatrix(result, TRUE)
  result@ownerPost <- ownerToMatrix(result, FALSE)

  ## Calibrate marginal costs from observed prices and demand params
  result@mcPre <- calcMC(result, TRUE)
  result@mcPost <- calcMC(result, FALSE)

  if (demand == "AIDS") {
    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result, ...)
  }


  ## Use observed prices as pre-merger equilibrium and only solve post-merger prices
  ## calcMC above already calibrated MCs using the supplied prices.
  result@pricePre <- prices
  if (!supply %in% c("auction", "bargaining2nd")) {
    result@pricePost <- calcPrices(result, FALSE, subset = subset, ...)
  } else {
    result@pricePost <- calcPrices(result, FALSE, ...)
  }

  if (any(grepl("logit", demand, ignore.case = TRUE), na.rm = TRUE)) {
    result@mktSize <- insideSize / sum(calcShares(result))
  } else if (any(grepl("ces", demand, ignore.case = TRUE), na.rm = TRUE)) {
    result@mktSize <- insideSize * (1 + result@slopes$alpha)
  }


  return(result)
}
