# VIRAD: Research Appendices and Supplementary Methods

Supplementary research documentation for **VIRAD**, a leakage-controlled multimodal framework for viral RNA motif discovery and Cas13d guide prioritization.

This repository preserves the detailed methods, algorithm descriptions, mathematical equations, supporting tables, computational environment, and supplementary figure legends from the main [VIRAD-2026 research repository](https://github.com/Ilovethesmiths1982/VIRAD-2026).

## Appendix navigation

| Document | Contents |
| --- | --- |
| [Detailed Methodology](DETAILED_METHODOLOGY.md) | Full research workflow, data acquisition, quality control, CD-HIT-EST cluster splitting, evolutionary conservation, 512-nt sliding windows, multimodal features, Viroformer architecture, training, decoding, motif evidence, guide design, ranking, and statistical methods |
| [Supplementary Methods and Algorithms](SUPPLEMENTARY_METHODS.md) | Sections S1–S20 with implementation-level details, mathematical equations, model architecture, leakage controls, Cas13d validation, panel freeze, and reproducibility checkpoints |
| [Appendices and Supplementary Tables](APPENDICES_AND_TABLES.md) | Appendices A–O: dataset and cluster counts, model tensors, training and validation results, decoder parameters, motif catalog, Cas13d guide evidence, frozen-panel test coverage, and supporting data index |
| [Computational Environment](COMPUTATIONAL_ENVIRONMENT.md) | Software versions, hardware, optimization configuration, random seeds, and reproducibility settings |
| [Supplementary Figure Legends](SUPPLEMENTARY_FIGURE_LEGENDS.md) | Captions and provenance notes for supplementary attention visualizations |

## Research-data archive

The source data are organized in the [VIRAD supporting-data Google Drive folder](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew), including motif catalogs, 300-guide development results, the frozen 18-guide panel, and untouched test evaluations. Availability to outside readers depends on Google Drive sharing settings.

## Reporting boundaries

- Motif discovery spans ten viral families; the Cas13d guide panel spans six families.
- Training, validation, held-out motif-model evaluation, and untouched frozen-guide-panel testing are distinct analyses.
- Predicted motifs and sequence-matched Cas13d guides are computational candidates, not experimentally demonstrated antiviral agents.
- The original CSV, FASTA, JSON, and figure artifacts are separate from these Markdown documents; the text alone does not replace those data files.

**Primary manuscript:** [RC27 VIRAD Article](https://github.com/Ilovethesmiths1982/VIRAD-2026/blob/main/manuscripts/RC27_VIRAD_ARTICLE_DRAFT.md).

**Source provenance:** These documents were copied from the main research repository on 2026-11-11. The original scientific source documents were retained without rewriting their findings.
