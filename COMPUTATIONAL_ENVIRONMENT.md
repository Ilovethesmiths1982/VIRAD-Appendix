# Appendix L. Computational Environment and Reproducibility

## L.1 Computational Environment

The VIRAD computational workflow was implemented primarily in Python 3.10 using PyTorch 2.0.1 with CUDA 11.8 support. The associated deep-learning environment included torchvision 0.15.2, torchaudio 2.0.2, Biopython 1.81, and torch-geometric 2.3.1. Model training, multimodal feature processing, inference, and downstream computational analyses were performed within controlled Conda environments so that the principal software dependencies remained consistent throughout the study.

Sequence redundancy reduction was performed using CD-HIT version 4.8.1. CD-HIT was developed for efficient clustering of highly similar biological sequences and was used in VIRAD before partition assignment so that closely related viral genomes remained within the same partition unit (Fu et al., 2012). Multiple-sequence alignment for evolutionary conservation analysis was performed using MAFFT version 7.520, which provides a widely used framework for scalable and accurate multiple-sequence alignment (Katoh & Standley, 2013). RNA-family annotation support was generated using Infernal version 1.1.5 together with Rfam covariance models. Rfam provides curated RNA-family models, while Infernal enables covariance-model searches that incorporate both sequence and RNA secondary-structure information (Kalvari et al., 2021; Nawrocki & Eddy, 2013).

Evolutionary modeling used the PHAST utilities phyloFit and phyloP. PHAST provides comparative-genomic tools for fitting phylogenetic models and calculating site-specific conservation or acceleration statistics (Hubisz et al., 2011). In VIRAD, phyloFit used a REV substitution model, and phyloP was applied in likelihood-ratio-test mode with CONACC scoring. The exact PHAST package release used during the completed analysis was not independently recorded and is therefore not assigned retrospectively.

RNA secondary-structure representations were generated using UFold, a deep-learning method for RNA secondary-structure prediction (Fu et al., 2022). VIRAD retained UFold-derived symmetric binary nucleotide-contact matrices as the RNA-structure modality. These matrices represented predicted contacts as discrete 0/1 values rather than calibrated base-pairing probabilities. Structure outputs were validated for expected matrix dimensions, finite values, binary semantics, matrix symmetry, and coordinate consistency before they were accepted for model input.

Host-context representations were generated using the ESM-2 `esm2_t33_650M_UR50D` protein language model. ESM-2 belongs to the evolutionary-scale protein language-model family described by Lin et al. (2023). The retained host representation had a dimensionality of 1,280 and was projected into VIRAD's shared 768-dimensional hidden representation before multimodal integration. These host embeddings were treated as global contextual representations and were not interpreted as direct nucleotide-level measurements of physical virus-host interaction.

For Cas13d candidate characterization, RNA structural and thermodynamic analyses were performed using ViennaRNA version 2.4.7. The ViennaRNA package provides computational tools for RNA secondary-structure prediction and thermodynamic analysis (Lorenz et al., 2011). RNAplfold was used to estimate target accessibility, RNAfold was used to evaluate spacer and complete-crRNA secondary structure, and RNAduplex was used to estimate guide-target interaction energies. Human-transcript similarity screening was performed using Bowtie2 against the GENCODE v49 human transcriptome. Bowtie2 is a fast gapped-read alignment method suitable for efficient sequence-similarity searches (Langmead & Salzberg, 2012). The exact Bowtie2 executable release used for the final screen was not independently recorded and is therefore not stated.

## L.2 Hardware and Training Runtime

The completed VIRAD training run was performed using an NVIDIA L40 GPU with approximately 46 GB of available GPU memory. A physical micro-batch size of four samples was combined with eight gradient-accumulation steps, resulting in an effective batch size of 32. Training used eight data-loader workers.

Most network operations were executed using BF16 mixed precision, while TensorFloat-32 operations and scaled dot-product attention were enabled where supported. Numerically sensitive sequence-to-structure and host cross-attention operations were locally evaluated in FP32 before their outputs were converted back to the working precision. These numerical settings were used to improve computational stability and efficiency without changing the model architecture, training data, optimization objective, or effective batch semantics.

The final training run continued for 51 epochs and reached 95,064 optimizer steps. Validation was performed after every epoch. The highest validation F1 was 0.8112 at epoch 50, while the final epoch produced a validation F1 of 0.8100. The complete epoch-by-epoch training and validation progression is reported separately in the training-performance appendix.

## L.3 Optimization Configuration

Optimization used AdamW with β1 = 0.9, β2 = 0.999, weight decay of 0.01, and global gradient-norm clipping at 1.0. AdamW decouples weight decay from the gradient-based parameter update and was used as the primary optimization algorithm in VIRAD (Loshchilov & Hutter, 2019). A learning-rate warmup period was followed by cosine annealing toward a minimum learning rate of 1 × 10^-6.

Training followed a four-stage curriculum in which biological modalities were introduced progressively. Stage 1 used sequence information only with an initial learning rate of 1 × 10^-4. Stage 2 introduced evolutionary conservation and RNA-structure information while retaining the same initial learning rate. Stage 3 added host-context information and used an initial learning rate of 5 × 10^-5. Stage 4 performed full multimodal and multitask optimization using an initial learning rate of 1 × 10^-5.

Dropout was set to 0.1 and stochastic depth to 0.10. Stochastic depth was introduced as a regularization strategy in which residual transformations can be randomly bypassed during training (Huang et al., 2016). An exponential-moving-average teacher with decay 0.999 was also maintained during training.

## L.4 Deterministic and Reproducible Analysis Settings

Deterministic procedures were used where random ordering or sampling could otherwise alter the analysis. Dataset partitioning used random seed 42 and assigned complete CD-HIT-EST clusters, rather than individual sequences or windows, to training, validation, and test partitions. This design reduced the risk that highly related viral genomes would be distributed across multiple partitions.

The final processed dataset contained 59,630 training windows from 1,724 clusters, 12,852 validation windows from 372 clusters, and 13,209 test windows from 375 clusters. No CD-HIT-EST cluster identifier was shared among the three partitions. Evolutionary conservation models were fitted using training genomes only. Validation and test genomes were introduced only after the training-defined alignment, evolutionary-distance model, neighbor-joining tree, and fitted substitution model had been established, thereby preventing held-out sequences from influencing the conservation feature background.

The interval-decoding procedure was also fixed before final downstream evaluation. Nucleotide probabilities were smoothed using a three-position window. High-confidence seeds required smoothed probabilities greater than 0.50, connected continuation positions required values greater than 0.40, and neighboring intervals separated by no more than 24 nucleotides were merged. Interval boundaries could subsequently expand through adjacent positions with smoothed probabilities greater than 0.20. Retained intervals were required to contain at least 20 nucleotides and to have a mean smoothed probability of at least 0.30. Touching or overlapping intervals were remerged after expansion.

## L.5 Statistical Reproducibility

Statistical analyses were performed using the biological grouping appropriate to the question being evaluated rather than treating overlapping nucleotide windows as independent observations. Sequence-, accession-, cluster-, or family-level grouping was retained according to the analysis.

Where bootstrap uncertainty was evaluated, 10,000 family-stratified paired bootstrap replicates were generated using random seed 42. Ninety-five percent percentile intervals were calculated using the 2.5th and 97.5th percentiles of the bootstrap replicate distribution.

For analyses involving multiple formal statistical comparisons, either the Benjamini-Hochberg false-discovery-rate procedure or the Holm sequential family-wise error procedure was used as appropriate. Benjamini and Hochberg (1995) introduced the widely used false-discovery-rate framework, whereas Holm (1979) described a sequentially rejective multiple-testing procedure. The correction method was defined with the corresponding analysis rather than selected retrospectively according to statistical significance.

---
