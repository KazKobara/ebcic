# Change Log

All notable changes, such as backward incompatibilities, will be documented in this file.

<!-- markdownlint-disable MD024 no-duplicate-heading -->
<!-- ## [Unreleased] -->

## [0.0.3]

### Added

- Approximation using the beta function.
- Params.exact_border() and Params.exact().

### Changed

- Max of x-axis in a graph to (`k_end * 10**log_n_end`) from (`3 * k_end * 10**log_n_end`).
- Made plot settings in interval_graph() in common as set_graph().
- For k=0, interval_graph() shows only upper intervals due to set_graph().
- Made k=0 and k=n be one-sided for 'wilson', 'wilson_cc' and 'beta_approx' in accordance with 'exact'.

## [0.0.2]

### Added

- Usage in MATLAB.
- Notes that 'confidence percentage' of k=0 or k=n is for one-sided.

## [0.0.1]

### Added

- Command line program using 'ebcic' package.
- 'ebcic' package is available with `pip install ebcic`

<!--
## Template
### Added
### Changed
### Deprecated
### Removed
### Fixed
### Security
-->
