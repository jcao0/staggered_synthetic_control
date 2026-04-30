# ---------- print ----------------------------------------------------------

#' Print an ssc Object
#'
#' Compact one-line summary of a \code{"ssc"} estimation result.
#'
#' @param x An object of class \code{"ssc"}, as returned by \code{\link{ssc}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.ssc <- function(x, ...) {
  cat("Staggered Synthetic Control\n")
  cat(sprintf("  N = %d units, T = %d pre-treatment, S = %d post-treatment\n",
              x$N, x$T, x$S))
  cat(sprintf("  K = %d treated (unit, period) pairs\n", x$K))
  cat(sprintf("  Design matrix min eigenvalue: %.4f\n", x$min_eigenvalue))

  cat(sprintf("\n  Overall ATT: %.4f", x$att_overall))
  if (!is.na(x$ci_lower_overall)) {
    cat(sprintf("  [%.4f, %.4f]  (%.0f%% CI)",
                x$ci_lower_overall, x$ci_upper_overall,
                100 * (1 - x$alpha)))
  }
  if (!is.na(x$p_value))
    cat(sprintf("  p = %.4f", x$p_value))
  cat("\n")
  invisible(x)
}


# ---------- summary --------------------------------------------------------

#' Summarise an ssc Object
#'
#' Prints a detailed summary of a \code{"ssc"} estimation result, including
#' design diagnostics, the overall ATT with confidence interval and p-value,
#' and a table of event-time ATT estimates.
#'
#' @param object An object of class \code{"ssc"}, as returned by
#'   \code{\link{ssc}}.
#' @param ... Currently unused.
#'
#' @return \code{object}, invisibly.
#'
#' @export
summary.ssc <- function(object, ...) {

  x <- object
  cat("Staggered Synthetic Control -- Summary\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("  Units (N):              %d\n",   x$N))
  cat(sprintf("  Pre-treatment (T):      %d\n",   x$T))
  cat(sprintf("  Post-treatment (S):     %d\n",   x$S))
  cat(sprintf("  Treated pairs (K):      %d\n",   x$K))
  cat(sprintf("  T / N ratio:            %.2f\n",  x$T / x$N))
  cat(sprintf("  Min eigenvalue:         %.6f\n", x$min_eigenvalue))
  cat(sprintf("  Significance level:     %.2f\n", x$alpha))

  has_ci <- !is.na(x$ci_lower_overall)

  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("  Overall ATT:            %.4f\n", x$att_overall))
  if (has_ci) {
    cat(sprintf("  %.0f%% CI:                [%.4f, %.4f]\n",
                100 * (1 - x$alpha), x$ci_lower_overall, x$ci_upper_overall))
    cat(sprintf("  p-value:                %.4f\n", x$p_value))
  } else {
    cat("  CI / p-value:           NA  (T < S)\n")
  }

  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("  Event-time ATT:\n\n")

  # build results table
  et <- seq_len(x$S) - 1L
  tab <- data.frame(
    event_time = et,
    att        = x$att_event,
    ci_lower   = x$ci_lower_event,
    ci_upper   = x$ci_upper_event
  )
  # only show first/last rows if S is large
  if (x$S > 20) {
    show <- rbind(utils::head(tab, 10), utils::tail(tab, 5))
    print(show, row.names = FALSE, digits = 4)
    cat(sprintf("  ... (%d rows total)\n", x$S))
  } else {
    print(tab, row.names = FALSE, digits = 4)
  }

  invisible(x)
}


# ---------- plot -----------------------------------------------------------

#' Plot Event-Time ATT from SSC Estimation
#'
#' @param x An object of class \code{"ssc"}.
#' @param main Title string.
#' @param xlab,ylab Axis labels.
#' @param ci Logical: draw the confidence band?  Default \code{TRUE} if
#'   inference is available.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object (invisibly) if \pkg{ggplot2} is available;
#'   otherwise a base-R plot is drawn and \code{NULL} is returned invisibly.
#'
#' @examples
#' set.seed(1)
#' N <- 10; Ttot <- 8
#' Y <- matrix(rnorm(N * Ttot), N, Ttot)
#' D <- matrix(0L, N, Ttot)
#' D[1:3, 5:Ttot] <- 1L   # units 1-3 treated from period 5
#' fit <- ssc(Y, D, S = 2, alpha = 0.05)
#' plot(fit)
#'
#' @export
plot.ssc <- function(x, main = "Event-time ATT (SSC)",
                     xlab = "Event time", ylab = "ATT estimate",
                     ci = !anyNA(x$ci_lower_event), ...) {

  et <- seq_len(x$S) - 1L   # event time starting at 0

  if (requireNamespace("ggplot2", quietly = TRUE)) {

    # workaround: avoid R CMD check NOTE for aes() column names
    event_time <- att <- ci_lower <- ci_upper <- NULL

    df <- data.frame(event_time = et, att = x$att_event)
    if (ci) {
      df$ci_lower <- x$ci_lower_event
      df$ci_upper <- x$ci_upper_event
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = event_time, y = att)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "gray40") +
      ggplot2::geom_line(ggplot2::aes(colour = "SSC"), linewidth = 0.9) +
      ggplot2::geom_point(ggplot2::aes(colour = "SSC"), size = 1.5)

    if (ci) {
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper, fill = "SSC"),
          alpha = 0.15)
    }

    p <- p +
      ggplot2::scale_colour_manual(name = NULL,
        values = c("SSC" = "#E8593C", "GSC" = "steelblue")) +
      ggplot2::scale_fill_manual(name = NULL,
        values = c("SSC" = "#E8593C", "GSC" = "steelblue"),
        guide = "none") +
      ggplot2::labs(title = main, x = xlab, y = ylab) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 14)
      )

    print(p)
    return(invisible(p))

  } else {

    # fallback to base R
    ylim <- range(c(x$att_event, x$ci_lower_event, x$ci_upper_event),
                  na.rm = TRUE)
    plot(et, x$att_event, type = "l", col = "#E8593C", lwd = 2,
         xlab = xlab, ylab = ylab, main = main, ylim = ylim)
    graphics::abline(h = 0, lty = 2, col = "gray40")
    if (ci && !anyNA(x$ci_lower_event)) {
      graphics::polygon(
        c(et, rev(et)),
        c(x$ci_lower_event, rev(x$ci_upper_event)),
        col = grDevices::adjustcolor("#E8593C", alpha.f = 0.15),
        border = NA
      )
      graphics::lines(et, x$att_event, col = "#E8593C", lwd = 2)
    }
    return(invisible(NULL))
  }
}
