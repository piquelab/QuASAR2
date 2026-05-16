## ============================================================================
##  testContrast.R - Wald test for a linear combination of coefficients.
## ============================================================================

#' Wald test for a general linear contrast \code{L * beta = 0}.
#'
#' @param object   a \code{quasar2CR} fit.
#' @param L        contrast vector or matrix.  Rows define contrasts (k = nrow);
#'                 column names must match \code{colnames(object$designs[[1]])}.
#'                 Use \code{\link{makeContrastsQuasar}} to build it from
#'                 string expressions.
#' @param ref      \code{"F"} (default) or \code{"chisq"}.  Use F at small
#'                 \code{n_g}; see \code{\link{testCoef}}.
#' @param df_method see \code{\link{testCoef}}.  Only used when \code{ref = "F"}.
#' @param df_extra numeric; only used when \code{df_method = "manual"}.
#' @return data frame with columns: identifier, statistic, df1, df2,
#'   p.value, padj.
#' @seealso \code{\link{testCoef}}, \code{\link{testLRT}},
#'   \code{\link{makeContrastsQuasar}}.
#' @export
testContrast <- function(object, L,
                         ref       = c("F", "chisq"),
                         df_method = c("residual", "moderated", "manual"),
                         df_extra  = 0) {
  stopifnot(inherits(object, "quasar2CR"))
  ref       <- match.arg(ref)
  df_method <- match.arg(df_method)
  sigma2    <- object$sigma2_psi
  if (is.null(sigma2)) sigma2 <- NA_real_
  if (is.null(dim(L))) L <- matrix(L, nrow = 1, dimnames = list(NULL, names(L)))

  rows <- lapply(seq_along(object$coefficients), function(k) {
    id <- names(object$coefficients)[k]
    b  <- object$coefficients[[k]]
    V  <- object$vcov[[k]]
    cn <- colnames(L);  if (is.null(cn)) cn <- names(b)
    miss <- setdiff(cn, names(b))
    if (length(miss))
      return(data.frame(identifier = id, statistic = NA, df1 = NA, df2 = NA,
                        p.value = NA))
    Lk  <- L[, names(b), drop = FALSE]
    est <- drop(Lk %*% b)
    VV  <- Lk %*% V %*% t(Lk)
    chi2 <- tryCatch(drop(crossprod(est, solve(VV, est))),
                     error = function(e) NA_real_)
    df1 <- nrow(Lk)
    n_g <- nrow(object$designs[[k]])
    p_g <- length(b)
    df2 <- switch(df_method,
                  residual  = n_g - p_g,
                  moderated = {
                    d0 <- if (is.finite(sigma2) && sigma2 > 0) 1 / sigma2 else 0
                    n_g - p_g + d0
                  },
                  manual    = n_g - p_g + df_extra)
    df2 <- max(1, df2)
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


#' Build a contrast matrix from string expressions of coefficient names.
#'
#' Inspired by \code{limma::makeContrasts} but pulls coefficient names from a
#' \code{quasar2CR} fit instead of a separate \code{levels} argument.
#'
#' @param contrasts character vector of expressions, e.g.
#'   \code{c("TreatmentB - TreatmentA", "TreatmentC - TreatmentA")}.
#'   Supports \code{+}, \code{-}, and \code{k*term} syntax.
#' @param object a \code{quasar2CR} fit.
#' @return matrix with one row per contrast expression, columns named by
#'   coefficient.
#' @examples
#' \dontrun{
#' L <- makeContrastsQuasar(c("trtB - trtA"), fit)
#' testContrast(fit, L)
#' }
#' @export
makeContrastsQuasar <- function(contrasts, object) {
  stopifnot(inherits(object, "quasar2CR"))
  cn <- colnames(object$designs[[1]])
  L  <- matrix(0, nrow = length(contrasts), ncol = length(cn),
               dimnames = list(contrasts, cn))
  for (i in seq_along(contrasts)) {
    expr <- contrasts[[i]]
    tokens <- regmatches(expr, gregexpr("[+-]?[^+-]+", expr))[[1]]
    tokens <- gsub("\\s", "", tokens)
    for (tk in tokens) {
      sign <- 1
      if (substr(tk, 1, 1) == "-") { sign <- -1; tk <- substring(tk, 2) }
      if (substr(tk, 1, 1) == "+") { tk <- substring(tk, 2) }
      mult <- 1
      if (grepl("\\*", tk)) {
        parts <- strsplit(tk, "\\*")[[1]]
        mult  <- as.numeric(parts[1])
        tk    <- parts[2]
      }
      if (!(tk %in% cn))
        stop("Contrast term '", tk, "' not in coefficient names: ",
             paste(cn, collapse = ", "))
      L[i, tk] <- L[i, tk] + sign * mult
    }
  }
  L
}
