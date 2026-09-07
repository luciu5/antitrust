#' @title Monopolistic-competition classes
#' @name MonCom-classes
#' @aliases MonComLogit-class MonComCES-class MonComBLP-class
#' @description `MonComLogit` and `MonComCES` contain the existing Logit and
#' CES demand state under a differentiated-product monopolistic-competition
#' conduct assumption.  Each product chooses its own price and does not
#' internalize strategic cross-product effects.  For `MonComCES`, the
#' aggregate CES index is held fixed in the perceived own-product derivative,
#' so the own elasticity used by the conduct FOC is `-gamma`; this is distinct
#' from the full share-adjusted elasticity returned by `elast()`.  For
#' `MonComBLP`, the derivative is integrated over the existing consumer draws
#' as `sum_r w_r * alpha_r * s_jr`.
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

#' @rdname MonCom-classes
#' @export
setClass(
  Class = "MonComBLP",
  contains = "LogitBLP"
)
