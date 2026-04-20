# Algorithm: NearRepeat_chunked()

This document explains the Knox test as implemented by Steenbeck's `NearRepeat()`
and how `NearRepeat_chunked()` performs the same computation with a different
memory strategy. The two implementations are mathematically equivalent.

---

## 1. The Knox test

The Knox test asks: are crime events that are spatially close also more likely
to be temporally close than would be expected by chance?

**Step 1 — Observed contingency table.**
For every ordered pair of events *(i, j)* with *i < j*, compute:
- Spatial distance: d_s(i, j) = Manhattan or Euclidean distance between coordinates
- Temporal distance: d_t(i, j) = |date_i − date_j| in days

Bin each distance into user-defined spatial and temporal bands (e.g. 0–1 m,
1–100 m, 101–200 m; same day, days 1–7, days 8–14). Count how many pairs fall
into each cell of the resulting *ns × nt* contingency table. This is the
**observed** table.

**Step 2 — Monte Carlo null distribution.**
Repeat `nrep` times:
1. Randomly permute the time labels across all events (spatial coordinates stay fixed).
2. Recount pairs by (spatial band, temporal band) under the shuffled times.

This builds a reference distribution of what the contingency table would look
like if space and time were independent.

**Step 3 — Knox ratio and p-value.**
For each cell (s, t):

```
Knox ratio  = observed[s, t] / mean(simulated[s, t, 1:nrep])
p-value     = (#{simulations where simulated ≥ observed} + 1) / (nrep + 1)
```

A Knox ratio substantially above 1.0 with a small p-value indicates that
crime pairs are over-represented in that spatio-temporal band — a near-repeat
signal. The conventional threshold used here is Knox ratio ≥ 1.2 and p < 0.05
(Ratcliffe's recommendation).

---

## 2. The standard approach: where `NearRepeat()` uses memory

Steenbeck's implementation computes distances using R's `stats::dist()`:

```r
# NearRepeat.R, lines 133-134
s_dist <- stats::dist(xy, method = method)   # spatial distances
t_dist <- stats::dist(mydf$time)             # temporal distances
```

`dist()` returns a *condensed distance matrix* — a vector of length *n(n−1)/2*
containing every pairwise distance. For *n* events, each value is a 64-bit
double (8 bytes):

```
RAM per matrix = n × (n − 1) / 2 × 8 bytes
```

Both matrices are alive simultaneously, so peak usage at this stage is:

```
~2 × n × (n − 1) / 2 × 8 bytes = n(n − 1) × 8 bytes
```

Inside the permutation loop (NearRepeat.R, line 155), a new temporal distance
matrix is created for each of the `nrep` shuffled time vectors:

```r
t_dist_perm <- stats::dist(sample(mydf$time))
```

This object is garbage-collected between iterations, so it does not accumulate,
but each iteration still peaks at one full matrix.

**Memory vs n:**

| n | Per matrix | Peak (observed phase) | Peak (per permutation) |
|---|---|---|---|
| 5,000 | 100 MB | 200 MB | 100 MB |
| 10,000 | 400 MB | 800 MB | 400 MB |
| 20,000 | 1.5 GB | 3 GB | 1.5 GB |
| 50,000 | 10 GB | 20 GB | 10 GB |
| 79,294 | 25 GB | **50 GB** | 25 GB |

---

## 3. The chunked approach: two-pass algorithm

`NearRepeat_chunked()` replaces the two `dist()` calls with explicit loops
over rows of the (never-materialised) distance matrix. At any point in the
computation, at most one row's worth of distances — a vector of length *n − i*
— exists in memory.

### Pass 1: build the observed table and collect within-range pairs

```
for i = 1 to n − 1:
    compute spatial distances from event i to all j > i
        (a vector of length n − i, not a full matrix)

    [optional] discard pairs where distance > max(sds)
        — valid when sds contains no Inf band

    bin remaining distances into spatial band indices with cut()

    for valid spatial pairs only:
        compute temporal distances: |time_i − time_j|
        bin into temporal band indices with cut()
        increment observed_vec[s_band + (t_band − 1) × ns]

    store (i, j, s_band) for each valid spatial pair
```

The spatial pre-filter is the key practical optimisation. In an urban crime
dataset with a maximum spatial band of 1001 m, the vast majority of event
pairs are further apart than this threshold and can be skipped entirely. For
Oslo residential burglary (~79k events spread across ~480 km²), this reduces
the number of pairs stored from ~3.1 billion to roughly 20–40 million —
a reduction of ~99%.

The pair list is built using pre-allocated list slots (one per row *i*) and
flattened with `unlist()` after the loop. This avoids the O(n²) memory
copying that would occur with repeated `c()` calls inside the loop.

### Pass 2: Monte Carlo permutations

```
for rep = 1 to nrep:
    perm_time = time_v[sample(1:n)]     ← shuffle time labels

    for each stored pair (pair_i, pair_j, s_band):
        t_dist = |perm_time[pair_i] − perm_time[pair_j]|
        bin into temporal band with cut()

    tabulate valid pairs into array_Knox[, , rep]
```

Because spatial distances are unchanged by the time permutation, the
(pair_i, pair_j, s_band) triples computed in Pass 1 are reused directly
across all `nrep` permutations. Only the temporal distance is recomputed,
and only for the within-range pairs.

Peak RAM during Pass 2 = storage for the pair list + one permuted time
vector + one array of pairwise temporal distances for the stored pairs.
For Oslo residential burglary: roughly 250–750 MB total.

---

## 4. Mathematical equivalence

The two implementations compute the same quantities.

**Observed counts.** Both implementations count the number of event pairs
*(i, j)* with *i < j* whose spatial distance falls in band *s* and temporal
distance falls in band *t*. The set of pairs considered, the distance metric,
and the `cut()` parameters (`right`, `include.lowest`, `dig.lab`) are
identical. The observed counts are therefore **exactly equal** — they are
fully deterministic with no stochastic element.

**Knox ratios and p-values.** Both use the same Monte Carlo procedure: shuffle
time labels, recount pairs under shuffled times, accumulate into `array_Knox`.
The p-value formula is identical:

```
p[s, t] = (#{rep : simulated[s, t, rep] ≥ observed[s, t]} + 1) / (nrep + 1)
```

The two functions will produce numerically different Knox ratios and p-values
in any given run because they use different random number sequences. This is
the same run-to-run variability that exists between two separate calls to
`NearRepeat()` — it is not a methodological difference. With nrep = 999 the
Monte Carlo error is small: the standard error of a p-value near 0.05 is
roughly ±0.007.

---

## 5. Temporal distance for Date objects

R's `stats::dist()` coerces Date vectors to their internal numeric
representation (integer days since 1970-01-01) before computing distances.
`NearRepeat_chunked()` performs the same coercion explicitly:

```r
time_v <- as.numeric(mydf$time)
```

The temporal distance between two events is then `|time_v[i] − time_v[j]|`,
which equals the number of calendar days between them. This is identical to
what `dist()` computes for Date input.

---

## 6. Band boundary convention

Both functions use `cut()` with `right = FALSE` and `include.lowest = FALSE`
(the defaults). This gives left-closed, right-open intervals:

```
sds = c(0, 1, 101, 201, ...)
→ bands: [0, 1),  [1, 101),  [101, 201), ...
```

With integer Manhattan distances, `[1, 101)` captures distances 1, 2, …, 100 m
exactly. The break at 101 rather than 100 is deliberate: a break at 100 would
produce `[1, 100)` which excludes distance = 100 m.

This convention follows Steenbeck's recommendation and avoids the off-by-one
labelling error present in the original Near Repeat Calculator (NRC) software,
where the label "a to b" corresponded to the computation interval [a−1, b).

---

## 7. Summary

| Property | NearRepeat() | NearRepeat_chunked() |
|---|---|---|
| Observed counts | Deterministic | **Identical** |
| Knox ratios / p-values | MC noise | MC noise (different seed) |
| Peak RAM | O(n²) | O(n + within_range_pairs) |
| Runtime | Fast (C-level dist) | Slower (R-level loop) |
| Feasible at n = 79k | No (crashes) | Yes |
| Output structure | `knox` class list | **Identical** |
