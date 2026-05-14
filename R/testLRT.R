## ============================================================================
##  testLRT.R - Likelihood ratio test against a reduced design.
## ============================================================================

#' Likelihood ratio test against a reduced design.
#'
#' For each SNP, refits \eqn{\hat\beta_g^{(r)}} on a reduced design,
#' **holding dispersion and eps fixed** at the full-model values, then
#' computes \eqn{\Lambda_g = 2(\ell_{full} - \ell_{red})}.
#'
#' Holding \eqn{\psi_g} fixed at the full-model MAP follows the
#' edgeR/DESeq2 convention and avoids spurious LRT inflation from
#' re-estimating dispersion under the null.
#'
#' By default the reference distribution is \code{F(df1, df2)} with
#' \code{df1 = p_full - p_red} and \code{df2} from \code{df_method} -
#' the small-n moderated analogue of \code{edgeR::glmQLFTest}.  Use
#' \code{ref = "chisq"} for the asymptotic chi-square reference; this
#' will be anti-conservative at small \code{n_g}.
#'
#' @param object   the full-design \code{quasar2CR} fit.
#' @param reduced  formula.  Must be a strict sub-set of
#'                 \code{object$formula}'s expanded columns.
#' @param ref      \code{"F"} (default) or \code{"chisq"}.
#' @param df_method see \code{\link{testCoef}}.
#' @param df_extra numeric; only when \code{df_method = "manual"}.
#' @return data frame with columns: identifier, statistic, df1, df2,
#'   p.value, padj.
#' @export
testLRT <- function(object, reduced,
                    ref        = c("F", "chisq"),
                    df_method  = c("residual", "moderated", "manual"),
                    df_extra   = 0) {
  stopifnot(inherits(object, "quasar2CR"))
  ref       <- match.arg(ref)
  df_method <- match.arg(df_method)
  sigma2    <- object$sigma2_psi
  if (is.null(sigma2)) sigma2 <- NA_real_
  dd <- object$data

  X_red_full <- stats::model.matrix(reduced, data = dd)

  rows <- lapply(seq_along(object$coefficients), function(k) {
    id    <- names(object$coefficients)[k]
    g     <- dd[dd$identifier == id, , drop = FALSE]
    X_red <- X_red_full[g$.row, , drop = FALSE]
    kp <- colnames(X_red) %in% "(Intercept)" |
          apply(X_red, 2, function(z) length(unique(z)) > 1)
    X_red <- X_red[, kp, drop = FALSE]

    psi_g <- object$dispersion$psi[k]
    fit_red <- tryCatch(
      .bb_fit_beta(g$R, g$A, X_red, g$offset, psi_g, object$eps),
      error = function(e) NULL
    )
    if (is.null(fit_red))
      return(data.frame(identifier = id, statistic = NA,
                        df1 = NA, df2 = NA, p.value = NA))
    ll_red  <- fit_red$value
    ll_full <- object$loglik[k]
    df_full <- length(object$coefficients[[k]])
    df_red  <- length(fit_red$beta)
    chi2 <- 2 * (ll_full - ll_red)
    chi2 <- pmax(chi2, 0)
    df1  <- df_full - df_red
    n_g  <- nrow(object$designs[[k]])
    df2 <- switch(df_method,
                  residual  = n_g - df_full,
                  moderated = {
                    d0 <- if (is.finite(sigma2) && sigma2 > 0) 1 / sigma2 else 0
                    n_g - df_full + d0
                  },
                  manual    = n_g - df_full + df_extra)
    df2 <- max(1, df2)
    if (df1 <= 0)
      return(data.frame(identifier = id, statistic = NA,
                        df1 = df1, df2 = df2, p.value = NA))
    if (ref == "F") {
      Fstat <- chi2 / df1
      p     <- stats::pf(Fstat, df1, df2, lower.tail = FALSE)
      data.frame(identifier = id, statistic = Fstat,
                 df1 = df1, df2 = df2, p.value = p, stringsAsFactors = FALSE)
    } else {
      p <- stats::pchisq(chi2, df1, lower.tail = FALSE)
      data.frame(identifier = id, statistic = chi2,
                 df1 = df1, df2 = NA, p.value = p, stringsAsFactors = FALSE)
    }
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$p.value, method = "BH")
  out
}
