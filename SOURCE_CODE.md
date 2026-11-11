# Source code excerpts and provenance

This directory contains a **selected snapshot** of implementation files copied without scientific modification from [VIRAD-2026](https://github.com/Ilovethesmiths1982/VIRAD-2026) on 2026-11-11.

## Included scripts

| Script | Purpose |
| --- | --- |
| [preprocess.py](scripts/preprocess.py) | Data preparation entry point |
| [train.py](scripts/train.py) | Training entry point |
| [validate.py](scripts/validate.py) | Validation and evaluation entry point |
| [infer.py](scripts/infer.py) | Model inference entry point |
| [validation_interval_grid.py](scripts/validation_interval_grid.py) | Interval-decoder validation analysis |
| [cas13d_full_validation.py](scripts/cas13d_full_validation.py) | Guide development, freeze integrity, and untouched-test evaluation |

## Selected model and data modules

| Module | Purpose |
| --- | --- |
| [Viroformer block](src/virad/models/viroformer.py) | Model architecture implementation |
| [VIRAD model](src/virad/models/virad.py) | Model assembly and prediction logic |
| [Partitioning](src/virad/data/splits.py) | Data splitting logic |
| [Sliding windows](src/virad/features/windows.py) | Window construction logic |

**Important:** This is a documentation-oriented **partial code snapshot**, not a complete standalone installable package. Scripts can import additional modules, expect configurations, external executables, model checkpoints, and input datasets that are not mirrored here. For runnable source, dependency installation, test suites, and environment files, use the [main repository](https://github.com/Ilovethesmiths1982/VIRAD-2026) and its [package metadata](https://github.com/Ilovethesmiths1982/VIRAD-2026/blob/main/pyproject.toml). No claim is made that these excerpts independently reproduce the reported results.

The included [MIT license](LICENSE) was copied from the original repository. External datasets, pretrained models, software, and third-party dependencies retain their respective licenses and attribution requirements.
