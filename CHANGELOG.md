# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] — 2026-07-26

### Added
- Initial release with complete RNA-seq differential expression workflow
- GEO data download script (`01_download_data.R`)
- Metadata preparation script (`02_prepare_metadata.R`)
- DESeq2 differential expression analysis script (`03_deseq2_analysis.R`)
- Visualization script (`04_visualization.R`)
- Functional enrichment analysis script (`05_functional_enrichment.R`)
- Shared utility functions (`utils.R`)
- Publication-quality figures (PCA, volcano, heatmap, MA plot, enrichment plots, pathview)
- GO enrichment (BP, MF, CC) with dotplots and barplots
- KEGG pathway enrichment with enrichment maps and cnetplots
- Pathway visualization via pathview (hsa05164 — Coronavirus disease)
- Scientific report (docs/report/) with introduction, methods, results, discussion, conclusion
- Presentation slides (docs/presentation/) with 6 chapters
- GitHub Pages website structure
- Issue templates, PR template, contribution guide, code of conduct, security policy
- MIT License
- CITATION.cff with DOI references
- GitHub Actions CI workflow
- badges for R, Bioconductor, DESeq2, license, commits, stars, issues, downloads

### Changed
- N/A (initial release)

### Fixed
- N/A (initial release)
