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

## APA 7 In-Text Citation Forms Used in This Appendix

| Source | Parenthetical citation | Narrative citation |
| --- | --- | --- |
| CD-HIT | (Fu et al., 2012) | Fu et al. (2012) |
| MAFFT | (Katoh & Standley, 2013) | Katoh and Standley (2013) |
| Rfam | (Kalvari et al., 2021) | Kalvari et al. (2021) |
| Infernal | (Nawrocki & Eddy, 2013) | Nawrocki and Eddy (2013) |
| PHAST / phyloP | (Hubisz et al., 2011) | Hubisz et al. (2011) |
| UFold | (Fu et al., 2022) | Fu et al. (2022) |
| ESM-2 | (Lin et al., 2023) | Lin et al. (2023) |
| ViennaRNA | (Lorenz et al., 2011) | Lorenz et al. (2011) |
| Bowtie2 | (Langmead & Salzberg, 2012) | Langmead and Salzberg (2012) |
| AdamW | (Loshchilov & Hutter, 2019) | Loshchilov and Hutter (2019) |
| Stochastic depth | (Huang et al., 2016) | Huang et al. (2016) |
| False-discovery-rate correction | (Benjamini & Hochberg, 1995) | Benjamini and Hochberg (1995) |
| Holm correction | (Holm, 1979) | Holm (1979) |

In APA 7 style, an ampersand (`&`) is used between two authors inside a parenthetical citation, whereas the word `and` is used when the authors are incorporated into the sentence as a narrative citation. For works with three or more authors, `et al.` is used from the first in-text citation.

---

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological), 57*(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: Accelerated for clustering the next-generation sequencing data. *Bioinformatics, 28*(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565

Fu, L., Cao, Y., Wu, J., Peng, Q., Nie, Q., & Xie, X. (2022). UFold: Fast and accurate RNA secondary structure prediction with deep learning. *Nucleic Acids Research, 50*(3), e14. https://doi.org/10.1093/nar/gkab1074

Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics, 6*, 65–70.

Huang, G., Sun, Y., Liu, Z., Sedra, D., & Weinberger, K. Q. (2016). Deep networks with stochastic depth. In B. Leibe, J. Matas, N. Sebe, & M. Welling (Eds.), *Computer vision—ECCV 2016* (pp. 646–661). Springer. https://doi.org/10.1007/978-3-319-46493-0_39

Hubisz, M. J., Pollard, K. S., & Siepel, A. (2011). PHAST and RPHAST: Phylogenetic analysis with space/time models. *Briefings in Bioinformatics, 12*(1), 41–51. https://doi.org/10.1093/bib/bbq072

Kalvari, I., Nawrocki, E. P., Ontiveros-Palacios, N., Argasinska, J., Lamkiewicz, K., Marz, M., Griffiths-Jones, S., Toffano-Nioche, C., Gautheret, D., Weinberg, Z., Rivas, E., Eddy, S. R., Finn, R. D., Bateman, A., & Petrov, A. I. (2021). Rfam 14: Expanded coverage of metagenomic, viral and microRNA families. *Nucleic Acids Research, 49*(D1), D192–D200. https://doi.org/10.1093/nar/gkaa1047

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. *Molecular Biology and Evolution, 30*(4), 772–780. https://doi.org/10.1093/molbev/mst010

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods, 9*(4), 357–359. https://doi.org/10.1038/nmeth.1923

Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., Smetanin, N., Verkuil, R., Kabeli, O., Shmueli, Y., dos Santos Costa, A., Fazel-Zarandi, M., Sercu, T., Candido, S., & Rives, A. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science, 379*(6637), 1123–1130. https://doi.org/10.1126/science.ade2574

Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology, 6*, 26. https://doi.org/10.1186/1748-7188-6-26

Loshchilov, I., & Hutter, F. (2019). Decoupled weight decay regularization. In *International Conference on Learning Representations*.

Nawrocki, E. P., & Eddy, S. R. (2013). Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics, 29*(22), 2933–2935. https://doi.org/10.1093/bioinformatics/btt509
