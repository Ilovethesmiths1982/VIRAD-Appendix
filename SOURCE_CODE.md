# Supplementary Source Code

This section documents selected VIRAD implementation files associated with the computational methodology.

## Workflow scripts

| File | Function |
| --- | --- |
| [Preprocessing](scripts/preprocess.py) | Sequence preprocessing entry point |
| [Training](scripts/train.py) | Model training entry point |
| [Validation](scripts/validate.py) | Model validation and evaluation |
| [Inference](scripts/infer.py) | Model inference |
| [Interval decoder analysis](scripts/validation_interval_grid.py) | Validation-based interval decoding |
| [Cas13d guide validation](scripts/cas13d_full_validation.py) | Guide development, frozen-panel selection, and untouched-test evaluation |

## Model and data modules

| File | Function |
| --- | --- |
| [Viroformer](src/virad/models/viroformer.py) | Transformer architecture components |
| [VIRAD model](src/virad/models/virad.py) | Model integration and prediction |
| [Dataset partitioning](src/virad/data/splits.py) | Dataset split implementation |
| [Sliding windows](src/virad/features/windows.py) | Window enumeration and representation |

## Reproducibility notes

These files constitute a **selected code excerpt**, not a complete installable or independently verified software package. Execution may require additional Python modules, configuration files, pretrained assets, external software, and input datasets.

Research input and evaluation files are indexed in the [VIRAD Supplementary Data archive](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew).

Code is distributed under the [MIT License](LICENSE). Third-party dependencies retain their original licenses.
