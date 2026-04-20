# generate_oslo_sim.R
# Creates the synthetic Oslo example dataset for the NearRepeatLite package.
#
# The real Oslo crime data (STRASAK police register) is not publicly available
# and is not distributed with this package. This script generates a fully
# synthetic stand-in that shares the same spatial and temporal extent.
#
# Spatial bounding box: UTM zone 32N (EPSG:25832), approximate Oslo extent
#   x: 258000 -- 275000 m (easting)
#   y: 6640000 -- 6665000 m (northing)
#
# Temporal range: 2015-01-01 to 2019-12-31 (five-year study period)
#
# To make the Knox test output informative as a demonstration, ~50 secondary
# events are injected within 200 m and 14 days of randomly chosen seed events.
# This produces mild near-repeat clustering without being based on real records.
#
# Run this script once from the package root to regenerate data/oslo_sim.rda
# and inst/extdata/oslo_sim.csv.
# -------------------------------------------------------------------------

set.seed(1519)

# Base events (random) -----------------------------------------------------
n_base   <- 450
x_base   <- runif(n_base, min = 258000, max = 275000)
y_base   <- runif(n_base, min = 6640000, max = 6665000)
all_days <- seq.Date(as.Date("2015-01-01"), as.Date("2019-12-31"), by = "day")
date_base <- sample(all_days, n_base, replace = TRUE)

# Near-repeat secondary events (mild clustering) ----------------------------
n_sec    <- 50
seed_idx <- sample(n_base, n_sec, replace = TRUE)

x_sec    <- x_base[seed_idx]    + runif(n_sec, -200, 200)
y_sec    <- y_base[seed_idx]    + runif(n_sec, -200, 200)
date_sec <- date_base[seed_idx] + sample(1:14, n_sec, replace = TRUE)

# Clamp secondary dates to study period
date_sec <- pmin(date_sec, as.Date("2019-12-31"))

# Combine and sort ----------------------------------------------------------
oslo_sim <- data.frame(
  x    = c(x_base, x_sec),
  y    = c(y_base, y_sec),
  date = c(date_base, date_sec)
)
oslo_sim <- oslo_sim[order(oslo_sim$date), ]
rownames(oslo_sim) <- NULL

# Save as package data (.rda) and raw CSV -----------------------------------
save(oslo_sim, file = "data/oslo_sim.rda", compress = "xz")
write.csv(oslo_sim, file = "inst/extdata/oslo_sim.csv", row.names = FALSE)

message("oslo_sim saved: ", nrow(oslo_sim), " events, ",
        min(oslo_sim$date), " to ", max(oslo_sim$date))
