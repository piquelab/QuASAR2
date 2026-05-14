## ============================================================================
##  methods-quasar2CR.R   -   S3 methods for class "quasar2CR".
## ============================================================================

#' Extract regression coefficients from a quasar2CR fit.
#' @param object a \code{quasar2CR} object.
#' @param id     optional character; if supplied, return only this identifier's
#'               coefficients.  Otherwise return the named list of all.
#' @param ...    ignored.
#' @return numeric vector (one SNP) or named list (all SNPs).
#' @export
coef.quasar2CR <- function(object, id = NULL, ...) {
  if (is.null(id)) return(object$coefficients)
  object$coefficients[[id]]
}

#' Extract Wald variance-covariance matrices from a quasar2CR fit.
#' @param object a \code{quasar2CR} object.
#' @param id     optional character; identifier of a single SNP.
#' @param ...    ignored.
#' @return matrix (one SNP) or named list (all SNPs).
#' @export
vcov.quasar2CR <- function(object, id = NULL, ...) {
  if (is.null(id)) return(object$vcov)
  object$vcov[[id]]
}

#' @export
print.quasar2CR <- function(x, ...) {
  cat("Quasar2 (Cox-Reid + EB) fit\n")
  cat("  Identifiers:    ", length(x$coefficients), "\n")
  cat("  Design:         ", deparse(x$formula), "\n")
  cat("  Coefficients:   ", paste(colnames(x$designs[[1]]), collapse = ", "), "\n")
  cat("  Shared eps:     ", signif(x$eps, 4), "\n")
  cat("  EB prior var psi:", signif(x$sigma2_psi, 4), "\n")
  cat("  Iterations:     ", x$n_iter, " (converged: ", x$converged, ")\n", sep = "")
  invisible(x)
}

#' Per-SNP tidy summary of a quasar2CR fit.
#'
#' Returns a long data frame with one row per (identifier, coefficient),
#' carrying point estimates, Wald standard errors, t-statistics and
#' p-values (using residual df = n_g - p_g).
#' For multi-coefficient testing use \code{\link{testContrast}}, and for
#' nested-model testing use \code{\link{testLRT}}.
#' @param object a \code{quasar2CR} fit.
#' @param ...    ignored.
#' @export
summary.quasar2CR <- function(object, ...) {
  rows <- lapply(seq_along(object$coefficients), function(k) {
    id  <- names(object$coefficients)[k]
    b   <- object$coefficients[[k]]
    V   <- object$vcov[[k]]
    se  <- sqrt(diag(V))
    n_g <- nrow(object$designs[[k]])
    p_g <- length(b)
    df  <- max(1, n_g - p_g)
    stat <- b / se
    p    <- 2 * stats::pt(-abs(stat), df = df)
    data.frame(identifier = id,
               term       = names(b),
               estimate   = unname(b),
               std.error  = unname(se),
               statistic  = unname(stat),
               df         = df,
               p.value    = unname(p),
               row.names  = NULL,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$p.value, method = "BH")
  out
}
