# panoramic 0.99.4

- Fixed local competition enrichment neighbor counts. `closepairs()` returns
  both orientations of a pair, and PANORAMIC now retains each pair once before
  adding it to the neighborhoods of its two endpoints.

- Added `stat = "local_comp_global_enrichment"` as an explicit alternative
  local-composition statistic. It uses a sample-wide target proportion
  (a random-label null conditional on observed locations), while the corrected
  edge-corrected CSR statistic remains the default.

- Fixed block-bootstrap handling: empty tiles no longer contribute artificial
  zero-valued block means, block-weighted estimates are invariant to point
  input order, and `boot = "block"` defaults to a 4-by-4 tile grid when tile
  counts are not supplied.

- Hardened input validation for sample identifiers, preparation parameters,
  and degenerate convex-hull windows.

- Build-system compatibility:
  - Switched the vignette to the `knitr::knitr` engine and added
    `markdown` to `Suggests` so `R CMD build` does not require Pandoc on
    Bioconductor builders.
- Documentation updates:
  - Updated `DESCRIPTION`, `README.md`, `inst/CITATION`, and the tutorial
    vignette to reflect the published *Bioinformatics* article
    (`10.1093/bioinformatics/btag546`) and current Bioconductor
    installation guidance.
  - Revised the package description to identify edge-corrected local
    competition enrichment as PANORAMIC's primary spatial statistic.
  - Removed generated help pages for non-exported implementation details so
    the installed manual exposes only the supported API.

# panoramic 0.99.3

- Major workflow/API restructuring:
  - Added `panoramic_analyze()` as the canonical end-to-end workflow.
  - `panoramic()` now delegates to `panoramic_analyze()` and returns a unified
    list with `prep`, `stats`, `pooled`, and `tables`.
  - Removed legacy compatibility arguments from the main wrappers
    (`meta_model`, `tau2`, `group_method`, `keep_boot`).
- Replaced legacy two-stage group-comparison flow with multilevel modeling:
  - Retired/unexported `panoramic_meta()` and `panoramic_compare_groups()`.
  - Added/exported `panoramic_meta_mv()` as the primary pooling API, including
    single-stage two-group contrasts (`beta_diff`, `se_diff`, `z_diff`,
    `p_diff`, `fdr_diff`) and optional group-specific heterogeneity estimates.
  - Contrast direction is now deterministic for grouped analyses.
- Added new result access helpers:
  - `panoramic_extract_contrast()` for extracting `spatialstats`, `meta`, or
    `contrast` tables.
  - Pre-flattened tables are stored in `metadata(se)$panoramic$tables`.
- Retired unused interpret/diagnostic helper prototypes
  (`panoramic_sample_deviation()`, `panoramic_interpret_meta_mv()`) from the
  public API to reduce surface area.
- Expanded statistical and plotting functionality:
  - Refactored `panoramic_spatialstats()` internals for robustness (single-sample
    support, clearer per-sample/pair failure handling, and optimized local
    composition enrichment caching/bootstrapping).
  - L-function variance is computed via a delta-method transform from
    corresponding K-function estimates.
  - Added `plot_representative_samples()` and updated `plot_volcano()`,
    `plot_forest()`, and network plotting helpers to align with multilevel
    contrast outputs.
- Documentation/submission updates:
  - Reworked vignette/tutorial to match the restructured workflow and current
    plotting APIs.
  - Added/updated man pages for newly exported functions and removed stale
    legacy docs.
  - Updated package metadata for Bioconductor review (`Depends: R (>= 4.6)`,
    CITATION parsing fixes, `.gitignore` final newline).
- Review-driven hardening and consistency fixes:
  - `panoramic_prepare()` now performs explicit schema validation
    (`design$sample`/`design$group`, duplicate sample IDs) and stronger
    coordinate/window checks after filtering.
  - Auto pair generation now uses observed marks (not unused factor levels),
    preventing inflated impossible pair sets under `min_presence`.
  - `panoramic_meta_mv()` now stops on ambiguous group labels after
    sanitization (`make.names`) to avoid silent group collisions.
  - `panoramic_spatialstats()` now degrades to per-sample NA outputs when local
    composition precompute fails, instead of aborting the full run.
  - `create_spatial_network()` now validates required rowData columns and
    documents/guards nonstandard `sig_operator = "gt"` behavior.
  - Plotting docs/messages now state the current two-group constraint
    (`plot_forest()`, `plot_representative_samples()`), and network man pages
    are synced to current function signatures.

# panoramic 0.99.0

- First Bioconductor submission of PANORAMIC.
- Added a full tutorial vignette showing how to prepare data, run pooled spatial analysis, and interpret results.
- Improved plotting functions for differential colocalization summaries.
- Added package citation metadata to make it easier to cite PANORAMIC and its associated preprint/release record.
