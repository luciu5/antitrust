#' @title Monopolistic-competition classes
#' @name MonCom-classes
#' @aliases MonComLogit-class MonComCES-class
#' @description `MonComLogit` and `MonComCES` contain the existing Logit and
#' CES demand state under a differentiated-product monopolistic-competition
#' conduct assumption.  Each product chooses its own price and does not
#' internalize strategic cross-product effects.
#' @section Extends:
#' `MonComLogit` extends `Logit`; `MonComCES` extends `CES`.
#' @include BertrandRUMClasses.R
NULL

#' @rdname MonCom-classes
#' @export
setClass(
  Class = "MonComLogit",
  contains = "Logit",
  prototype = prototype(
    control.slopes = list(reltol = .Machine$double.eps^0.25)
  )
)

#' @rdname MonCom-classes
#' @export
setClass(
  Class = "MonComCES",
  contains = "CES",
  prototype = prototype(
    control.slopes = list(reltol = .Machine$double.eps^0.25)
  )
)
