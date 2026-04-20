# NearRepeatLite 0.1.0

* Initial release.
* Implements the Knox test for near-repeat spatio-temporal clustering using
  row-by-row distance computation instead of a full pairwise distance matrix.
* Reduces peak RAM from O(n²) to O(n + within-range pairs). For a dataset of
  ~79,000 events, peak usage drops from approximately 50 GB to under 1 GB.
* Output is fully compatible with the `knox` class and `plot.knox()` from
  Steenbeck's NearRepeat package.
* Developed for a near-repeat analysis of crime in Oslo, Norway (2015--2019).
