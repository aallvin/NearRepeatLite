# verify_equivalence.R
# Optional verification script: confirms that NearRepeatLite() and
# NearRepeat() produce identical observed counts for the same data and bands.
#
# Requirements:
#   - NearRepeatLite installed (remotes::install_github("aallvin/NearRepeatLite"))
#   - NearRepeat installed (remotes::install_github("wsteenbeek/NearRepeat"))
#
# Knox ratios and p-values will differ between the two functions due to Monte
# Carlo noise — this is expected and not a sign of error. Only the observed
# counts (which are fully deterministic) must be identical.
# -------------------------------------------------------------------------

library(NearRepeatLite)
library(NearRepeat)

# Load synthetic example data -----------------------------------------------
data(oslo_sim)

# Band definitions -----------------------------------------------------------
sds <- c(0, 1, 101, 201, 301, 401, 501)
tds <- c(0, 1, 8, 15, 22, 29)

# Run NearRepeatLite() -------------------------------------------------------
set.seed(1519)
result_lite <- NearRepeatLite(
  x    = oslo_sim$x,
  y    = oslo_sim$y,
  time = oslo_sim$date,
  sds  = sds,
  tds  = tds,
  nrep = 99
)

# Run NearRepeat() -----------------------------------------------------------
set.seed(1519)
result_std <- NearRepeat(
  x    = oslo_sim$x,
  y    = oslo_sim$y,
  time = oslo_sim$date,
  sds  = sds,
  tds  = tds,
  nrep = 99
)

# Compare observed counts ----------------------------------------------------
obs_lite <- as.matrix(result_lite$observed)
obs_std  <- as.matrix(result_std$observed)
diff_mat <- obs_lite - obs_std

cat("\n--- Observed counts: NearRepeatLite ---\n")
print(obs_lite)

cat("\n--- Observed counts: NearRepeat (Steenbeck) ---\n")
print(obs_std)

cat("\n--- Difference (should be all zeros) ---\n")
print(diff_mat)

if (all(diff_mat == 0)) {
  cat("\nVERIFICATION PASSED: observed counts are identical.\n")
} else {
  cat("\nVERIFICATION FAILED: observed counts differ.\n")
  cat("Non-zero cells:\n")
  print(which(diff_mat != 0, arr.ind = TRUE))
}
