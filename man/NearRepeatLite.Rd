\name{NearRepeatLite}
\alias{NearRepeatLite}
\title{Memory-Efficient Knox Test for Near-Repeat Spatio-Temporal Clustering}
\usage{
NearRepeatLite(x, y, time, sds, tds,
               s_include.lowest = FALSE, s_right = FALSE,
               t_include.lowest = FALSE, t_right = FALSE,
               method = "manhattan", nrep = 999,
               saveSimulations = FALSE, future.seed = TRUE, ...)
}
\arguments{
  \item{x}{Numeric vector of x-coordinates (e.g. easting in metres).}
  \item{y}{Numeric vector of y-coordinates (e.g. northing in metres).
    Must be the same length as \code{x} and \code{time}.}
  \item{time}{Vector of event dates or times. Coerced to numeric (integer
    days since 1970-01-01 for \code{Date} objects). Must be the same length
    as \code{x} and \code{y}.}
  \item{sds}{Numeric vector of spatial distance band breakpoints (at least
    two values). Intervals are left-closed, right-open \eqn{[a, b)} by
    default (\code{s_right = FALSE}).}
  \item{tds}{Numeric vector of temporal distance band breakpoints (at least
    two values). Same interval convention as \code{sds}.}
  \item{s_include.lowest}{Logical; whether the leftmost spatial interval
    should also be closed on the right. Default \code{FALSE}.}
  \item{s_right}{Logical; if \code{TRUE}, spatial intervals are right-closed
    \eqn{(a, b]}; if \code{FALSE} (default), left-closed \eqn{[a, b)}.}
  \item{t_include.lowest}{Logical; whether the leftmost temporal interval
    should also be closed on the right. Default \code{FALSE}.}
  \item{t_right}{Logical; if \code{TRUE}, temporal intervals are right-closed;
    if \code{FALSE} (default), left-closed.}
  \item{method}{Distance metric: \code{"manhattan"} (default) or
    \code{"euclidean"}.}
  \item{nrep}{Integer; number of Monte Carlo permutations. Must be \eqn{\ge 2}.
    Use 99 for exploratory runs and 999 for publication-ready results.
    Default \code{999}.}
  \item{saveSimulations}{Logical; if \code{TRUE}, the full simulation array
    is included in the output as \code{array_Knox}. Default \code{FALSE}.}
  \item{future.seed}{Accepted for compatibility with \code{NearRepeat()};
    not used. Reproducibility is controlled by \code{set.seed()} before
    the call.}
  \item{...}{Additional arguments accepted for compatibility; not used.}
}
\value{
A list of class \code{"knox"} with the following elements:
\describe{
  \item{observed}{An \eqn{n_s \times n_t} contingency table of observed
    crime-pair counts per spatial-temporal band combination.}
  \item{knox_ratio}{An \eqn{n_s \times n_t} matrix of Knox ratios:
    observed count divided by the mean of the simulated counts.}
  \item{knox_ratio_median}{An \eqn{n_s \times n_t} matrix of Knox ratios
    using the median of the simulated counts.}
  \item{pvalues}{An \eqn{n_s \times n_t} matrix of Monte Carlo p-values,
    computed as \eqn{(\#\{sim \geq obs\} + 1) / (nrep + 1)}.}
  \item{array_Knox}{An \eqn{n_s \times n_t \times nrep} array of simulated
    counts. Only present when \code{saveSimulations = TRUE}.}
}
The output is fully compatible with \code{plot.knox()} from the
\pkg{NearRepeat} package.
}
\description{
A drop-in replacement for \code{NearRepeat()} from Steenbeck's
\pkg{NearRepeat} package that uses row-by-row distance computation instead
of a full pairwise distance matrix. Implements the same Knox test for
spatio-temporal clustering but reduces peak RAM from \eqn{O(n^2)} to
\eqn{O(n + \text{within-range pairs})}.

For a dataset of ~79,000 events, peak RAM drops from approximately 50 GB
to under 1 GB. Results are identical in observed counts; Knox ratios and
p-values differ only due to Monte Carlo noise, as expected.
}
\examples{
data(oslo_sim)

sds <- c(0, 1, 101, 201, 301, 401, 501)
tds <- c(0, 1, 8, 15, 22, 29)

set.seed(1519)
result <- NearRepeatLite(
  x    = oslo_sim$x,
  y    = oslo_sim$y,
  time = oslo_sim$date,
  sds  = sds,
  tds  = tds,
  nrep = 99
)

result$observed
result$knox_ratio
result$pvalues
}
\seealso{
\code{vignette("NearRepeatLite")} for a worked example.
}
