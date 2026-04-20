#' Memory-Efficient Knox Test for Near-Repeat Spatio-Temporal Clustering
#'
#' @description
#' A drop-in replacement for \code{NearRepeat()} from Steenbeck's
#' \pkg{NearRepeat} package that uses row-by-row distance computation instead
#' of a full pairwise distance matrix. Implements the same Knox test for
#' spatio-temporal clustering but reduces peak RAM from \eqn{O(n^2)} to
#' \eqn{O(n + \text{within-range pairs})}.
#'
#' For a dataset of ~79,000 events, peak RAM drops from approximately 50 GB
#' (infeasible on standard hardware) to under 1 GB. Results are identical in
#' observed counts; Knox ratios and p-values differ only due to Monte Carlo
#' noise, as expected.
#'
#' See \code{vignette("NearRepeatLite")} for a worked example and
#' \code{algorithm.md} in the package source for a technical walkthrough of
#' both algorithms and their mathematical equivalence.
#'
#' @param x Numeric vector of x-coordinates (e.g. easting in metres).
#' @param y Numeric vector of y-coordinates (e.g. northing in metres).
#'   Must be the same length as \code{x} and \code{time}.
#' @param time Vector of event dates or times. Coerced to numeric (integer
#'   days since 1970-01-01 for \code{Date} objects). Must be the same length
#'   as \code{x} and \code{y}.
#' @param sds Numeric vector of spatial distance band breakpoints (at least
#'   two values). Intervals are left-closed, right-open \eqn{[a, b)} by
#'   default (\code{s_right = FALSE}). Example:
#'   \code{c(0, 1, 101, 201, 301, 401, 501, 601, 701, 801, 901, 1001)}.
#' @param tds Numeric vector of temporal distance band breakpoints (at least
#'   two values). Same interval convention as \code{sds}. Example:
#'   \code{c(0, 1, 8, 15, 22, 29)} for same-day and four weekly bands.
#' @param s_include.lowest Logical; whether the leftmost spatial interval
#'   should be closed on the right as well. Default \code{FALSE}.
#' @param s_right Logical; if \code{TRUE}, spatial intervals are
#'   right-closed (a, b]; if \code{FALSE} (default), left-closed \eqn{[a, b)}.
#' @param t_include.lowest Logical; whether the leftmost temporal interval
#'   should be closed on the right as well. Default \code{FALSE}.
#' @param t_right Logical; if \code{TRUE}, temporal intervals are
#'   right-closed; if \code{FALSE} (default), left-closed.
#' @param method Distance metric: \code{"manhattan"} (default) or
#'   \code{"euclidean"}.
#' @param nrep Integer; number of Monte Carlo permutations. Must be \eqn{\ge 2}.
#'   Use 99 for exploratory runs and 999 for publication-ready results.
#'   Default \code{999}.
#' @param saveSimulations Logical; if \code{TRUE}, the full
#'   \eqn{n_s \times n_t \times} \code{nrep} simulation array is included in
#'   the output as \code{array_Knox}. Default \code{FALSE}.
#' @param future.seed Accepted for compatibility with \code{NearRepeat()};
#'   not used. Reproducibility is controlled by \code{set.seed()} before
#'   the call.
#' @param ... Additional arguments accepted for compatibility; not used.
#'
#' @return A list of class \code{"knox"} with the following elements:
#' \describe{
#'   \item{observed}{An \eqn{n_s \times n_t} contingency table of observed
#'     crime-pair counts per spatial-temporal band combination.}
#'   \item{knox_ratio}{An \eqn{n_s \times n_t} matrix of Knox ratios:
#'     observed count divided by the mean of the simulated counts.}
#'   \item{knox_ratio_median}{An \eqn{n_s \times n_t} matrix of Knox ratios
#'     using the median of the simulated counts.}
#'   \item{pvalues}{An \eqn{n_s \times n_t} matrix of Monte Carlo p-values,
#'     computed as \eqn{(\#\{sim \geq obs\} + 1) / (nrep + 1)}.}
#'   \item{array_Knox}{An \eqn{n_s \times n_t \times nrep} array of simulated
#'     counts. Only present when \code{saveSimulations = TRUE}.}
#' }
#' The output is fully compatible with \code{plot.knox()} from the
#' \pkg{NearRepeat} package.
#'
#' @seealso
#' \code{\link[NearRepeat]{NearRepeat}} for the standard implementation;
#' \code{vignette("NearRepeatLite")} for a worked example.
#'
#' @examples
#' # Load the synthetic Oslo example dataset
#' data(oslo_sim)
#'
#' # Define spatial bands: repeat (0 m) + 100 m intervals to 500 m
#' sds <- c(0, 1, 101, 201, 301, 401, 501)
#'
#' # Define temporal bands: same day + four weekly bands
#' tds <- c(0, 1, 8, 15, 22, 29)
#'
#' # Run the Knox test (use nrep = 99 for a quick example)
#' set.seed(1519)
#' result <- NearRepeatLite(
#'   x    = oslo_sim$x,
#'   y    = oslo_sim$y,
#'   time = oslo_sim$date,
#'   sds  = sds,
#'   tds  = tds,
#'   nrep = 99
#' )
#'
#' # Inspect results
#' result$observed
#' result$knox_ratio
#' result$pvalues
#'
#' @importFrom progressr progressor with_progress
#' @export
NearRepeatLite <- function(x, y, time, sds, tds,
                           s_include.lowest = FALSE, s_right = FALSE,
                           t_include.lowest = FALSE, t_right = FALSE,
                           method = "manhattan", nrep = 999,
                           saveSimulations = FALSE, future.seed = TRUE, ...) {

  # Input validation ---------------------------------------------------------
  if (nrep < 2) stop("nrep must be >= 2")
  if (length(sds) < 2) stop("sds must have at least 2 elements")
  if (length(tds) < 2) stop("tds must have at least 2 elements")

  # Build and clean data frame -----------------------------------------------
  df <- data.frame(x = x, y = y, time = as.numeric(time))
  df <- df[stats::complete.cases(df), ]
  x_v    <- df$x
  y_v    <- df$y
  time_v <- df$time
  n  <- nrow(df)
  ns <- length(sds) - 1
  nt <- length(tds) - 1

  # Pre-compute band labels via cut() on empty vector (avoids repeated calls)
  s_labels <- levels(cut(numeric(0), breaks = sds,
                         right = s_right, include.lowest = s_include.lowest))
  t_labels <- levels(cut(numeric(0), breaks = tds,
                         right = t_right, include.lowest = t_include.lowest))

  # Spatial pre-filter threshold (skip pairs beyond max finite band boundary)
  max_sds <- if (is.finite(max(sds))) max(sds) else Inf

  # Pass 1: observed counts and within-range pair collection -----------------
  observed_vec <- integer(ns * nt)

  # Pre-allocate one list slot per row to avoid quadratic memory from c()
  pair_s_list <- vector("list", n - 1L)
  pair_i_list <- vector("list", n - 1L)
  pair_j_list <- vector("list", n - 1L)

  with_progress({
    p <- progressor(steps = n - 1L)

    for (i in seq_len(n - 1L)) {
      j_idx <- seq.int(i + 1L, n)

      # Spatial distances for row i
      if (method == "manhattan") {
        s_dist <- abs(x_v[i] - x_v[j_idx]) + abs(y_v[i] - y_v[j_idx])
      } else {
        s_dist <- sqrt((x_v[i] - x_v[j_idx])^2 + (y_v[i] - y_v[j_idx])^2)
      }

      # Spatial pre-filter: discard pairs beyond the largest finite band
      if (is.finite(max_sds)) {
        keep <- s_dist < max_sds
        if (!any(keep)) {
          p()
          next
        }
        j_idx  <- j_idx[keep]
        s_dist <- s_dist[keep]
      }

      # Bin spatial distances
      s_bins <- as.integer(cut(s_dist, breaks = sds,
                               right = s_right,
                               include.lowest = s_include.lowest))

      valid_s <- !is.na(s_bins)
      if (!any(valid_s)) {
        p()
        next
      }
      j_valid  <- j_idx[valid_s]
      s_valid  <- s_bins[valid_s]

      # Temporal distances and binning for valid spatial pairs
      t_dist <- abs(time_v[i] - time_v[j_valid])
      t_bins <- as.integer(cut(t_dist, breaks = tds,
                               right = t_right,
                               include.lowest = t_include.lowest))

      valid_t <- !is.na(t_bins)
      if (any(valid_t)) {
        # Linear index into ns x nt matrix, then tabulate
        lin_idx <- s_valid[valid_t] + (t_bins[valid_t] - 1L) * ns
        observed_vec <- observed_vec +
          tabulate(lin_idx, nbins = ns * nt)
      }

      # Store within-range spatial pairs for permutation reuse
      pair_s_list[[i]] <- s_valid
      pair_i_list[[i]] <- rep.int(i, length(j_valid))
      pair_j_list[[i]] <- j_valid

      p()
    }
  })

  # Flatten pair lists
  pair_s_band <- unlist(pair_s_list, use.names = FALSE)
  pair_i_idx  <- unlist(pair_i_list, use.names = FALSE)
  pair_j_idx  <- unlist(pair_j_list, use.names = FALSE)

  # Build observed table
  observed <- matrix(observed_vec, nrow = ns, ncol = nt,
                     dimnames = list(s_labels, t_labels))
  observed <- as.table(observed)

  # Pass 2: Monte Carlo permutations -----------------------------------------
  n_pairs      <- length(pair_s_band)
  array_Knox   <- array(0L, dim = c(ns, nt, nrep))

  with_progress({
    p2 <- progressor(steps = nrep)

    for (rep in seq_len(nrep)) {
      perm_time <- time_v[sample.int(n)]
      t_perm    <- abs(perm_time[pair_i_idx] - perm_time[pair_j_idx])
      t_bins_p  <- as.integer(cut(t_perm, breaks = tds,
                                  right = t_right,
                                  include.lowest = t_include.lowest))
      valid_p   <- !is.na(t_bins_p)
      if (any(valid_p)) {
        lin_p <- pair_s_band[valid_p] + (t_bins_p[valid_p] - 1L) * ns
        array_Knox[, , rep] <- array_Knox[, , rep] +
          tabulate(lin_p, nbins = ns * nt)
      }
      p2()
    }
  })

  # Knox ratios and p-values -------------------------------------------------
  mean_sim   <- apply(array_Knox, 1:2, mean)
  median_sim <- apply(array_Knox, 1:2, stats::median)

  knox_ratio        <- observed / mean_sim
  knox_ratio_median <- observed / median_sim

  pvalues <- matrix(NA_real_, nrow = ns, ncol = nt,
                    dimnames = list(s_labels, t_labels))
  for (i in seq_len(ns)) {
    for (j in seq_len(nt)) {
      pvalues[i, j] <- (sum(array_Knox[i, j, ] >= observed[i, j]) + 1L) /
        (nrep + 1L)
    }
  }

  # Assemble output ----------------------------------------------------------
  out <- list(
    observed          = observed,
    knox_ratio        = knox_ratio,
    knox_ratio_median = knox_ratio_median,
    pvalues           = pvalues
  )
  if (saveSimulations) out$array_Knox <- array_Knox
  class(out) <- "knox"
  out
}
