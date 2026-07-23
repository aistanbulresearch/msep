# Changelog

All notable changes to `msep` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [1.2.0]

### Changed
- **Ferroptosis gene set standardized to the analysis set** (`msep.pathways.FERROPTOSIS`).
  Five genes were swapped so the package reproduces the published figures:
  added `HSPB1`, `IREB2`, `NFS1`, `PROM2`, `STEAP3`; removed `ALOX15B`, `CARS1`,
  `CISD2`, `PEBP1`, `SLC3A2`. This changes ferroptosis pathway CV values (e.g.
  chordoma CSC ferroptosis CV 5.387 → 5.61); pin `msep==1.1.0` if you need the
  previous set.
- **`IMMUNE_EVASION` de-duplicated**: `PVRL2` removed as a withdrawn HGNC alias of
  `NECTIN2`; the set now resolves to 28 distinct genes (was declared 29).
- **`coordination.pathway_cv` now reports the number of genes actually averaged**
  (positive mean) rather than the number merely found in the matrix. No CV value
  changes; only the returned `n_genes` count.

## [1.1.0]

- PyPI release; bundled synthetic demo dataset and Colab quickstart.

## [1.0.1]

- Metadata fixes.

## [1.0.0]

- Initial release: per-cell Shannon entropy, pathway-level CV, bootstrap CIs,
  pseudo-perturbation, and XBP1 consolidation analysis.
