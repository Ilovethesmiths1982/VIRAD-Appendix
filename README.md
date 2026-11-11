# VIRAD
### Supplementary Methods, Data, and Source Code

**Researcher:** Zhyrex G. Sablay  
**Institution:** Tarlac National High School

**VIRAD** is a leakage-controlled multimodal computational framework for identifying conserved viral RNA motifs and prioritizing candidate Cas13d guides. This repository presents the study's extended methodology, mathematical formulations, supplementary tables, computational documentation, and selected implementation files in a structured research appendix.

## Repository contents

| Section | Description |
| --- | --- |
| [Detailed Methodology](DETAILED_METHODOLOGY.md) | Data acquisition, sequence quality control, cluster-disjoint partitioning, evolutionary conservation, multimodal features, Viroformer architecture, optimization, inference, motif characterization, and Cas13d guide prioritization |
| [Supplementary Methods and Algorithms](SUPPLEMENTARY_METHODS.md) | Extended protocols, mathematical equations, architecture details, validation procedures, and reproducibility checkpoints |
| [Supplementary Tables and Appendices](APPENDICES_AND_TABLES.md) | Dataset composition, training and evaluation summaries, motif catalogs, guide screening, frozen-panel selection, and untouched-test results |
| [Computational Environment](COMPUTATIONAL_ENVIRONMENT.md) | Hardware, software, optimization settings, and reproducibility information |
| [Source Code](SOURCE_CODE.md) | Index of selected implementation files and their functions |
| [Reproducibility Protocol](REPRODUCIBILITY_PROTOCOL.md) | Required artifacts, example command interfaces, data-firewall checks, and independent verification checklist |

## Research workflow

1. **Data preparation:** viral sequence acquisition, quality filtering, redundancy reduction, and leakage-controlled partitioning.
2. **Feature construction:** sliding-window representation of sequence, evolutionary conservation, predicted RNA structure, and host-associated context.
3. **Motif discovery:** multimodal learning, validation-derived interval decoding, and genome-scale motif catalog construction.
4. **Candidate guide development:** Cas13d guide enumeration, sequence and structural assessment, host-transcript similarity checks, and ranking analyses.
5. **Frozen-panel evaluation:** prespecified guide selection followed by evaluation against an untouched test partition.

## Supplementary datasets

The accompanying research data are available through the following archive:

**[VIRAD Supplementary Data — Google Drive](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew)**

The archive is organized into four collections:

| Collection | Contents |
| --- | --- |
| Motif catalog | Identified motif intervals, annotation comparisons, family-level summaries, and unmatched candidates |
| Guide development | Candidate guides, development validation, matched controls, and ranking sensitivity |
| Frozen guide panel | Selected 18-guide panel and freeze documentation |
| Untouched test | Frozen-panel test evaluation, coverage tables, and summary statistics |

Access to these files depends on the Google Drive folder's sharing settings.

## Source code

Selected research scripts and model modules are included in the [source-code index](SOURCE_CODE.md), covering preprocessing, training, validation, inference, interval decoding, Cas13d guide assessment, model architecture, partitioning, and window construction.

The included code is a **partial research snapshot**, not a complete standalone software distribution. See the [reproducibility protocol](REPRODUCIBILITY_PROTOCOL.md) for documented command interfaces, required input manifests, and current limitations. Additional dependencies, input data, model assets, and configurations may be required for execution.

## Scientific scope

The reported motif and guide findings are based on computational analyses. Annotation-unmatched motifs are not experimentally validated novel motifs. Predicted guide coverage does not demonstrate antiviral efficacy or establish biological safety. Development-set analyses and untouched-test results are reported separately to preserve their respective evaluation boundaries.

## License

Source code is provided under the [MIT License](LICENSE). Third-party software, datasets, and pretrained models remain subject to their respective licenses and usage conditions.

---

*VIRAD Supplementary Research Archive · 2026*
