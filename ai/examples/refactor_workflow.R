## New calibration/simulation workflow on the refactor branch.

prices <- c(2, 2.2, 2.5)
shares <- c(.35, .25, .20)
margins <- c(.40, .35, .30)
ownerPre <- c("A", "B", "C")
ownerPost <- c("A", "A", "C")

## Existing one-call workflow remains available.
old_result <- logit(
    prices = prices, shares = shares, margins = margins,
    ownerPre = ownerPre, ownerPost = ownerPost
)

## Calibrate once, then simulate one or more counterfactuals.
fit <- calibrate(
    demand = "logit", conduct = "cournot",
    prices = prices, shares = shares, margins = margins,
    ownerPre = ownerPre
)
result_ab <- simulate(fit, ownerPost = ownerPost)
result_ac <- simulate(fit, ownerPost = c("A", "B", "B"))

## Supplied structural parameters use the same simulation boundary.
specified_fit <- specify(
    demand = "logit", conduct = "bertrand",
    prices = prices,
    parameters = list(alpha = -1.2, meanval = c(.5, .3, .1)),
    ownerPre = ownerPre,
    shares = shares, margins = margins, insideSize = 100
)
specified_result <- simulate(specified_fit, ownerPost = ownerPost)

## ALM is an explicit model-specific calibration variant.  Its legacy
## constructor estimates the unknown outside share from additional margins.
alm_fit <- calibrate(
    demand = "logit", conduct = "bertrand", variant = "alm",
    prices = prices, shares = shares / sum(shares), margins = margins,
    ownerPre = ownerPre
)
alm_result <- simulate(alm_fit, ownerPost = ownerPost)

## Cost efficiencies and product exits are scenario inputs.
efficiency_result <- simulate(
    fit, ownerPost = ownerPost,
    mcDelta = c(-.05, 0, 0), subset = c(TRUE, TRUE, FALSE)
)
