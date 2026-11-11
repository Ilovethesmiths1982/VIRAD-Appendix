# VIRAD | Supplementary Research Archive

**Methods · Algorithms · Equations · Supplementary Tables · Source Code**

This repository accompanies the **VIRAD** research manuscript and provides extended methodological documentation, supplementary results tables, computational provenance, and selected source-code files for a leakage-controlled multimodal framework for **viral RNA motif discovery and Cas13d guide prioritization**.

**Main research repository:** [VIRAD-2026](https://github.com/Ilovethesmiths1982/VIRAD-2026)  
**Manuscript:** [RC27 VIRAD Article Draft](https://github.com/Ilovethesmiths1982/VIRAD-2026/blob/main/manuscripts/RC27_VIRAD_ARTICLE_DRAFT.md)  
**Supporting data archive:** [VIRAD Google Drive](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew)

## About this archive

VIRAD combines sequence, evolutionary conservation, predicted RNA secondary structure, and host-associated computational features to identify candidate viral RNA motifs. Downstream analyses prioritize Cas13d guides and evaluate a frozen guide panel against an untouched sequence partition. The archive is intended to make the computational workflow, model assumptions, validation boundaries, and supplementary evidence easier to inspect.

## Documentation

| Resource | Scope |
| --- | --- |
| **[Detailed Methodology](DETAILED_METHODOLOGY.md)** | End-to-end study design, sequence acquisition, filtering, leakage-controlled partitions, feature construction, architecture, optimization, inference, guide ranking, and statistical procedures |
| **[Supplementary Methods](SUPPLEMENTARY_METHODS.md)** | Sections S1–S20 with mathematical equations, implementation-level explanations, Viroformer information flow, validation protocols, and audit checkpoints |
| **[Appendices and Tables](APPENDICES_AND_TABLES.md)** | Appendices A–O, including dataset composition, model configuration, motif catalog summaries, guide evidence, panel freeze, and untouched-test results |
| **[Computational Environment](COMPUTATIONAL_ENVIRONMENT.md)** | Computational resources, software environment, reproducibility settings, and implementation provenance |
| **[Supplementary Figure Legends](SUPPLEMENTARY_FIGURE_LEGENDS.md)** | Supplementary figure descriptions and provenance notes |
| **[Source Code Index](SOURCE_CODE.md)** | Selected original source files, their roles, and the limits of this partial code snapshot |

## Selected implementation files

| Research component | Source |
| --- | --- |
| Sequence preprocessing | [scripts/preprocess.py](scripts/preprocess.py) |
| Model training | [scripts/train.py](scripts/train.py) |
| Model validation | [scripts/validate.py](scripts/validate.py) |
| Inference | [scripts/infer.py](scripts/infer.py) |
| Validation-derived interval decoding | [scripts/validation_interval_grid.py](scripts/validation_interval_grid.py) |
| Cas13d guide development, freeze, and test | [scripts/cas13d_full_validation.py](scripts/cas13d_full_validation.py) |
| Viroformer architecture | [src/virad/models/viroformer.py](src/virad/models/viroformer.py) |
| VIRAD model | [src/virad/models/virad.py](src/virad/models/virad.py) |
| Data partitioning | [src/virad/data/splits.py](src/virad/data/splits.py) |
| Sliding-window features | [src/virad/features/windows.py](src/virad/features/windows.py) |

The included code is **not a standalone distribution**. For all package modules, dependencies, configuration files, and tests, refer to the [main implementation repository](https://github.com/Ilovethesmiths1982/VIRAD-2026).

## Supporting data

The [research-data folder](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew) organizes four supporting datasets:

1. **Motif catalog:** identified motifs, annotation overlap, family summaries, and unmatched candidates.
2. **300-guide development:** development validation, matched controls, ranking sensitivity, and guide FASTA.
3. **Frozen 18-guide panel:** selected guide table and freeze manifest.
4. **Untouched test:** frozen-panel evaluation, coverage tables, and test summary.

The data files are linked rather than duplicated here. Public availability depends on the Google Drive sharing permissions.

## Interpretation and limitations

VIRAD's reported findings are **computational**. Annotation-unmatched motifs are not experimentally confirmed novel motifs. Sequence-level guide matches do not establish cellular activity, antiviral efficacy, or clinical safety. Guide development and frozen-panel testing must remain separate when interpreting results. See the supplementary methods and tables for analysis-specific definitions, population sizes, and thresholds.

## Reproducibility and provenance

Documentation and selected implementation files were copied from [VIRAD-2026](https://github.com/Ilovethesmiths1982/VIRAD-2026) on **11 November 2026**. The supplementary archive preserves the scientific source text rather than silently revising reported results. Files in this repository do not constitute a complete independent reproduction package, and inclusion of a method does not imply experimental verification.

## Citation and license

When referencing the research, cite the accompanying manuscript and identify the repository version or commit used. No DOI is asserted for this archive.

Original code is distributed under the copied [MIT License](LICENSE). Third-party software, data, and pretrained models may have separate licenses and terms.

---

**Repository purpose:** transparent supplementary documentation for scientific review, research communication, and computational reproducibility.
