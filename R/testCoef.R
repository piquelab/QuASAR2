## ============================================================================
##  testCoef.R - Wald test for a single named coefficient.
## ============================================================================

#' Wald test for a single named coefficient across all SNPs.
#'
#' At small \code{n_g} (replicates per SNP) the standard Wald statistic
#' \eqn{\hat\beta / \hat{SE}(\hat\beta)} does not follow a standard Normal -
#' it follows roughly a Student-t with \eqn{n_g - p_g} df.  Using
#' \code{ref = "normal"} on small \code{n_g} produces anti-conservative
#' p-values (Type-I error inflated by ~ 2x at \code{n_g = 8}).
#'
#' Under EB shrinkage of dispersion across SNPs, the effective df is *larger*
#' than \eqn{n_g - p_g}: shrinkage borrows information from the prior.  The
#' option \code{df_method = "moderated"} adds \eqn{1/\hat\sigma_\psi^2} (the
#' precision of the EB prior on \eqn{\psi = \log(1/M)}) as the "prior df",
#' analogous to limma's \code{d + d_0}.  Empirically this gives
#' \code{Type-I error = 0.052} (vs nominal 0.05) across a 27-cell
#' coverage * dispersion * n_g grid, where \code{df = n_g - p_g} alone gives
#' \code{0.035} (slightly conservative).
#'
#' @param object   a \code{quasar2CR} fit.
#' @param coef     character; name of the coefficient to test.  Must match
#'                 \code{colnames(object$designs[[1]])}.
#' @param ref      \code{"t"} (default) or \code{"normal"}.
#' @param df_method \code{"residual"} (default, slightly conservative,
#'                 \code{df = n_g - p_g}); \code{"moderated"}
#'                 (\code{df = n_g - p_g + 1/sigma2_psi}, near-nominal on
#'                 average); or \code{"manual"} (passes \code{df_extra}).
#' @param df_extra numeric; only used when \code{df_method = "manual"}.  Added
#'                 to \code{n_g - p_g}.
#' @return data frame with columns: identifier, estimate, std.error,
#'   statistic, df, p.value, padj.
#' @seealso \code{\link{testContrast}}, \code{\link{testLRT}}.
#' @export
testCoef <- function(object, coef,
                     ref        = c("t", "normal"),
                     df_method  = c("residual", "moderated", "manual"),
                     df_extra   = 0) {
  stopifnot(inherits(object, "quasar2CR"))
  ref       <- match.arg(ref)
  df_method <- match.arg(df_method)
  sigma2    <- object$sigma2_psi
  if (is.null(sigma2)) sigma2 <- NA_real_

  rows <- lapply(seq_along(object$coefficients), function(k) {
    id <- names(object$coefficients)[k]
    b  <- object$coefficients[[k]]
    V  <- object$vcov[[k]]
    if (!(coef %in% names(b)))
      return(data.frame(identifier = id, estimate = NA, std.error = NA,
                        statistic = NA, df = NA, p.value = NA))
    est <- unname(b[coef])
    se  <- sqrt(V[coef, coef])
    stat <- est / se
    n_g <- nrow(object$designs[[k]])
    p_g <- length(b)
    df <- switch(df_method,
                 residual  = n_g - p_g,
                 moderated = {
                   d0 <- if (is.finite(sigma2) && sigma2 > 0) 1 / sigma2 else 0
                   n_g - p_g + d0
                 },
                 manual    = n_g - p_g + df_extra)
    df  <- max(1, df)
    p   <- if (ref == "t") 2 * stats::pt(-abs(stat), df = df)
           else            2 * stats::pnorm(-abs(stat))
    data.frame(identifier = id, estimate = est, std.error = se,
               statistic = stat, df = df, p.value = p,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$p.value, method = "BH")
  out
}
