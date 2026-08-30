# Synthetic, non-confidential markets used by the formal QA suite.
#
# Keep these values deliberately small and asymmetric.  They are fixtures, not
# empirical data, and are also used by the reproducibility runner when it
# records checksums.

qa_fixture_market <- function() {
    list(
        prices = c(2.00, 2.20, 2.50),
        shares = c(0.35, 0.25, 0.20),
        margins = c(0.40, 0.35, 0.30),
        ownerPre = c("A", "B", "C"),
        ownerPost = c("A", "A", "C"),
        mcDelta = c(0, 0, 0),
        subset = c(TRUE, TRUE, TRUE),
        labels = c("p1", "p2", "p3")
    )
}

qa_fixture_diversions <- function() {
    matrix(c(
        -1, 0.55, 0.45,
        0.50, -1, 0.50,
        0.40, 0.60, -1
    ), nrow = 3, byrow = TRUE)
}

qa_fixture_nests <- function() c("N1", "N1", "N2")

qa_fixture_blp_parameters <- function() {
    list(
        alpha = -1,
        sigma = 0.10,
        sigmaNest = 1,
        nDraws = 24,
        # Supplying draws makes the stochastic fixture deterministic without
        # changing the caller's RNG stream.
        consDraws = seq(-1, 1, length.out = 24)
    )
}

qa_fixture_linear_parameters <- function() {
    list(
        slopes = matrix(c(
            -1.20, 0.10, 0.08,
             0.09, -1.10, 0.12,
             0.07, 0.11, -1.00
        ), nrow = 3, byrow = TRUE),
        intercepts = c(3.0, 3.1, 3.2)
    )
}

qa_fixture_aids_parameters <- function() {
    list(
        slopes = matrix(c(
            -0.13, 0.08, 0.05,
             0.08, -0.14, 0.06,
             0.05, 0.06, -0.11
        ), nrow = 3, byrow = TRUE),
        intercepts = c(0.40, 0.35, 0.25),
        mktElast = -1.2
    )
}
