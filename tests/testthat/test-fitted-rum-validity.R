.fitted_rum_validity_wrapper <- function(fit, target) {
  model <- fit@model
  values <- lapply(slotNames(model), function(name) slot(model, name))
  names(values) <- slotNames(model)
  values$parmsStart <- c(-1, .2)
  if (identical(target, "LogitCapALM")) {
    values$capacitiesPre <- values$capacitiesPost <- rep(Inf, length(model@shares))
  }
  suppressWarnings(do.call(new, c(list(Class = target), values)))
}

test_that("fitted ALM wrappers retain known normalization and demand", {
  cases <- list(
    list(demand = "logit", conduct = "bertrand", target = "LogitALM"),
    list(demand = "logit", conduct = "cournot", target = "LogitCournotALM"),
    list(demand = "ces", conduct = "bertrand", target = "CESALM"),
    list(demand = "logit", conduct = "bertrand", target = "LogitCapALM"),
    list(demand = "logit", conduct = "auction2nd", target = "Auction2ndLogitALM")
  )
  for (case in cases) {
    parameters <- if (case$demand == "ces") {
      list(gamma = 3, alpha = .4, meanval = c(1, .8, 1.2))
    } else list(alpha = -1.5, meanval = c(1, .8, 1.2))
    source <- suppressWarnings(specify(
      case$demand, case$conduct, prices = c(2, 2.2, 2.5),
      parameters = parameters, ownerPre = c("A", "B", "C"),
      insideSize = 100
    ))
    expect_false(isTRUE(all.equal(source@model@shareInside, 1)))
    target <- .fitted_rum_validity_wrapper(source, case$target)
    expect_true(suppressWarnings(validObject(target)), info = case$target)
    expect_identical(target@shareInside, source@model@shareInside)
    expect_identical(target@slopes, source@model@slopes)
    expect_equal(calcShares(target), calcShares(source@model), tolerance = 1e-12)
    expect_equal(elast(target), elast(source@model), tolerance = 1e-12)

    ## Without fitted parameters, these remain calibration inputs and must
    ## satisfy the historical ALM inside-share marker.
    target@slopes <- list()
    expect_error(suppressWarnings(validObject(target)), "sum of 'shares' must equal 1")
  }
})

test_that("fitted single-product ALM state needs no identifying margin pair", {
  source <- suppressWarnings(specify(
    "logit", "bertrand", prices = 2,
    parameters = list(alpha = -1.5, meanval = 1),
    ownerPre = "A", insideSize = 100
  ))
  expect_false(isTRUE(all.equal(source@model@shareInside, 1)))
  target <- .fitted_rum_validity_wrapper(source, "LogitALM")
  expect_true(suppressWarnings(validObject(target)))
  expect_equal(calcShares(target), calcShares(source@model), tolerance = 1e-12)
  target@slopes <- list()
  expect_error(suppressWarnings(validObject(target)), "At least 2 elements")
})
