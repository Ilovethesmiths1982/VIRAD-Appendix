# VIRAD Reproducibility Protocol

This protocol distinguishes **documented methods**, **available code**, and **artifacts still required for an independent reproduction**. Commands below describe the interfaces of the archived scripts; they are **not claimed to have been executed or validated in this repository**.

## 1. Reproducibility levels

| Level | Current status | Requirements |
| --- | --- | --- |
| Methods inspection | Available | Detailed methodology, equations, algorithm descriptions, and supplementary tables |
| Selected implementation inspection | Available | Script and model excerpts, configuration files, environment specification |
| Data/results audit | Partially available | Supporting CSV, JSON, and FASTA files in the [data archive](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew); original split and accession manifests must be checked |
| Independent end-to-end rerun | Not demonstrated | Complete source dependency tree, exact input datasets, all external tool versions, frozen manifests, checkpoints, and executable validation |
| Experimental validation | Not performed | Biological guide activity, viral suppression, and host safety experiments |

## 2. Required inputs and provenance

An independent reproduction requires the following records before workflow execution:

- Accession/version identifiers, source database, retrieval date, raw viral FASTA/GFF3, and any RNA annotations.
- Sequence QC and CD-HIT-EST clustering output with cluster IDs and family assignments.
- Train/validation/test accession manifest and cluster-disjointness audit.
- Training-only alignment and conservation-model inputs, including tool versions and parameters.
- UFold and host-context feature inputs and their provenance.
- Exact training configuration, random seed, checkpoints, and validation-selected decoder parameters.
- Motif catalog and Cas13d development inputs, including frozen-panel manifest and untouched-test population.

The [supplementary data folder](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew) contains motif, guide-development, freeze, and test collections. **The archive has not been verified to contain every upstream input or a complete executable reproduction package.**

## 3. Software and configuration

The archived configuration and environment files are:

- [configs/significant-130.yaml](configs/significant-130.yaml)
- [configs/paper.yaml](configs/paper.yaml)
- [environment.yml](environment.yml)
- [pyproject.toml](pyproject.toml)
- [Computational environment](COMPUTATIONAL_ENVIRONMENT.md)

These files were preserved as supplied. Independent verification requires comparison with the frozen run before any parameter can be treated as the definitive production value. Applicable operating-system, Python, CUDA/driver, PyTorch, CD-HIT-EST, MAFFT, PHAST, UFold, Bowtie2, and ViennaRNA versions should be recorded.

**Caution:** The repository includes selected modules, not all dependencies imported by the scripts. Installation or execution may fail until the remaining implementation files and external programs are supplied.

## 4. Example command interfaces

The following commands illustrate the existing CLI arguments. Angle-bracket placeholders represent verified local file paths. They are **templates**, not a claim of successful execution.

### Preprocessing

```bash
python scripts/preprocess.py \
  --fasta <viral_sequences.fasta> \
  --out <processed_output_directory> \
  --config configs/significant-130.yaml \
  --annotations <annotations.gff3>
```

Production preprocessing requires consistent annotation coordinate conventions and the appropriate structure and conservation inputs. Reproduction of structure-dependent analyses requires structure prediction to remain enabled.

### Training

```bash
python scripts/train.py \
  --data <processed_data_directory> \
  --config configs/significant-130.yaml \
  --checkpoint-dir <checkpoint_directory> \
  --seed 42
```

A seed alone does not guarantee bitwise reproducibility across hardware, CUDA versions, or nondeterministic kernels. The actual run configuration and device must be recorded.

### Inference

```bash
python scripts/infer.py \
  --fasta <input_sequences.fasta> \
  --checkpoint <frozen_checkpoint> \
  --config configs/paper.yaml \
  --out-dir <inference_output_directory>
```

Feature inputs and model weights must correspond to the same trained model and data conventions.

### Cas13d development, freeze, and untouched test

```bash
python scripts/cas13d_full_validation.py develop \
  --guides <development_guides.fasta> \
  --viral-fasta <viral_population.fasta> \
  --manifest <partition_manifest.csv> \
  --outdir <development_output_directory>
```

```bash
python scripts/cas13d_full_validation.py freeze \
  --development-csv <development_validation.csv> \
  --panel-size-per-family 3 \
  --outdir <frozen_panel_output_directory>
```

```bash
python scripts/cas13d_full_validation.py test \
  --frozen-panel <frozen_panel.csv> \
  --freeze-json <panel_freeze.json> \
  --viral-fasta <viral_population.fasta> \
  --manifest <partition_manifest.csv> \
  --panel-mismatch-threshold 2 \
  --outdir <test_output_directory>
```

The actual development analysis may require additional optional inputs for family alignments, motif catalogs, human-transcriptome screening, guide structure, and activity evidence. The minimal command does not reproduce every optional analysis.

## 5. Independent validation checklist

- [ ] Cluster disjointness across training, validation, and test is independently recomputed and verified.
- [ ] Sequence identifiers and window coordinates are checked against the original manifest.
- [ ] Conservation fitting is verified to use training records only.
- [ ] Model and decoder selection are verified against the permitted development partitions.
- [ ] Reported and reproduced metrics are compared using identical evaluation units and thresholds.
- [ ] The frozen guide-panel file hash is verified before untouched-test evaluation.
- [ ] Test-sequence exclusion from guide generation, development ranking, and threshold selection is verified.
- [ ] Per-family and overall test coverage are recalculated from raw sequence matches.
- [ ] Deviations, unavailable artifacts, and unreproduced results are documented without substituting approximations.

## 6. Reported analysis checkpoints

These are **reported results to verify**, not automated pass/fail tests: 59,630 training windows, 12,852 validation windows, 13,209 test windows; 1,424 cataloged motifs; 300 development guide candidates; 18 frozen guides; and 106 of 243 untouched test accessions matched by the panel at a threshold of at most two mismatches.

These counts belong to different analytical stages and must not be treated as interchangeable sample sizes. Family-specific denominators, definitions, and limitations are documented in the supplementary tables.

## 7. Known limitations

The current public appendix does not establish that all raw data, model checkpoints, full package modules, intermediate artifacts, and external tool installations are available. The archived documentation may describe optional or historical analyses; an independent rerun should verify which branches were active in the frozen experiment. The results are computational and do not establish experimental antiviral efficacy or host safety.
