# RC27 VIRAD Appendices and Supplementary Tables

**Companion to:** [RC27 VIRAD Article](RC27_VIRAD_ARTICLE_DRAFT.md)

These appendices organize the verified numerical tables and supporting information from the original research paper and the repository's existing `APPENDICES.md`. They are designed to accompany the RC27 manuscript, not to introduce additional experiments. Detailed algorithms and equations appear in [RC27 Detailed Methods](RC27_DETAILED_METHODS.md); the [supplementary attention figure captions](RC27_SUPPLEMENTARY_FIGURE_LEGENDS.md) are maintained separately.

**Scientific reporting distinction:** The validation partition was used for development and model evaluation. The separate 243-accession assessment tested the frozen 18-guide panel. Their metrics must never be described as the same test. The ten-family motif-discovery analysis and six-family guide-panel analysis likewise use different scopes.

## Appendix contents

| Appendix | Supporting material |
|---|---|
| A | Frozen partition counts and leakage control |
| B | Sequence collection and RNA annotation |
| C | Training-only evolutionary conservation |
| D | Synchronized multimodal input tensors |
| E | Six-block Viroformer model architecture |
| F | Training trajectory and **validation** interval performance |
| G | Validation-frozen interval decoder |
| H | Genome-scale motif inference |
| I | Conservation, structure, and motif evidence |
| J | Cas13d guide design and development screens |
| K | Independent frozen-panel accession coverage |
| L | Software, computational environment, reproducibility |
| M | Motif catalog and prioritized unmatched subset |
| N | Guide-level evidence and supplementary data inventory |

The original figure exports, complete 1,424-row motif catalog, complete 300-guide development table, and full frozen 18-guide files remain separate data artifacts. A textual appendix does not substitute for those files.

---

## Appendix A. Dataset Composition and Leakage-Controlled Partitioning

Viral genomes from ten target families were filtered to remove sequences shorter than 500 nt or containing more than 10% ambiguous bases. Local windows with more than 20% ambiguous bases or gaps were excluded. Sequences were clustered with CD-HIT-EST at 95% nucleotide identity and 90% bidirectional coverage, and complete clusters,not individual windows,were assigned to training, validation, or test partitions to prevent leakage (Fu et al., 2012).

**Table A1**  
*Frozen processed dataset composition by partition*

| Partition | Windows | CD-HIT-EST clusters |
| --- | ---: | ---: |
| Training | 59,630 | 1,724 |
| Validation | 12,852 | 372 |
| Test (held-out) | 13,209 | 375 |
| **Total** | **85,691** | **2,471** |

No cluster identifier was shared across partitions.

---

## Appendix B. Viral Acquisition and RNA Annotation

Viral genomes were retrieved from NCBI Entrez nuccore, prioritizing RefSeq records where available (Sayers, 2022). RNA-relevant supervision combined NCBI GFF3 annotations with Rfam covariance-model evidence identified using Infernal (Kalvari et al., 2021; Nawrocki & Eddy, 2013). Motif presence was encoded as a nucleotide-level binary target, while motif identity used an 87-class vocabulary where valid labels were available. Internal intervals used zero-based, half-open coordinates `[start, end)`.

---

## Appendix C. Training-Only Evolutionary Conservation Modeling

Training genomes were aligned within each viral family using MAFFT 7.520 (Katoh & Standley, 2013). Canonical nucleotide mismatches were converted to Jukes-Cantor distances, followed by deterministic neighbor joining (Jukes & Cantor, 1969; Saitou & Nei, 1987). The training alignment and tree were analyzed with PHAST phyloFit using a REV model and phyloP using LRT-CONACC scoring (Hubisz et al., 2011). Validation and test sequences were projected only after the training-derived evolutionary background had been frozen; no refitting was performed.

---

## Appendix D. Multimodal Feature Construction

VIRAD used synchronized 512-nt windows with four biological modalities: nucleotide sequence, phyloP conservation, UFold-derived RNA structure, and ESM-2 host context. UFold outputs were retained as symmetric binary contact matrices rather than calibrated probabilities (Fu et al., 2022). Host context was generated using `esm2_t33_650M_UR50D` and projected from 1,280 to 768 dimensions (Lin et al., 2023).

**Table D1**  
*Model-ready multimodal input representations*

| Component | Representation |
| --- | --- |
| Viral sequence | `[B,5,L]` for A/U/G/C/N |
| Conservation | `[B,L]` |
| UFold contacts | `[B,L,L]` |
| Host context | `[B,1280]` before projection |
| Shared hidden state | `[B,L,768]` |

---

## Appendix E. VIRAD and Viroformer Architecture

The production model contained 68,819,416 trainable parameters and six Viroformer blocks with hidden width 768, 12 self-attention heads, feed-forward width 3,072, dropout 0.1, and stochastic depth 0.10. Sequence-to-structure cross-attention used six heads. Each block combined self-attention with learned conservation and structure biases, host cross-attention, residual normalization, and a 768→3,072→768 feed-forward network. A dense pair decoder, symmetry enforcement, one structural-refinement pass, and a learned structure gate preceded the nucleotide-level prediction heads. The attention architecture follows the Transformer framework of Vaswani et al. (2017). **In simple terms:** the six blocks repeatedly compare nucleotide positions while incorporating evolutionary conservation and predicted RNA pairing into their attention scores. This architecture does not prove that highly attended positions are functionally important. For the exact equations and seven-part implementation-level explanation, see [Detailed Methods, Section S19](RC27_DETAILED_METHODS.md) and the source-paper crosswalk in Section S20.

---

**Dataset and evaluation distinction.** Table A1 counts **CD-HIT-EST clusters** and overlapping **512-nt windows**, not necessarily unique genome accessions. Validation contains **372 clusters and 12,852 windows**, whereas held-out model testing contains **375 clusters and 13,209 windows**. The source's interval evaluation separately describes **372 validation sequences**. Equality of that sequence count and the cluster count does not establish a one-sequence-per-cluster relationship. The reported family-level interval F1 values are validation results and cannot be relabeled held-out test results. The **2,497 processed viral segments** used for the genome-scale motif catalog represent a different analysis stage and should not be equated with the **2,471 clusters**.

## Appendix F. Training and Validation Interval Evaluation

Training used a four-stage curriculum: sequence only; sequence plus conservation and structure; addition of host context; and full multimodal optimization. AdamW was used with β1=0.9, β2=0.999, weight decay 0.01, gradient clipping at 1.0, warmup, and cosine decay (Loshchilov & Hutter, 2019). The effective batch size was 32. The final run ended at epoch 51 and optimizer step 95,064.

**Table F1**  
*Key training milestones*

| Epoch | Training loss | Validation loss | Validation F1 |
| ---: | ---: | ---: | ---: |
| 1 | 0.8921 | 0.9412 | 0.5120 |
| 10 | 0.6761 | 0.7352 | 0.6410 |
| 20 | 0.5001 | 0.5603 | 0.7440 |
| 30 | 0.3755 | 0.4378 | 0.7980 |
| 40 | 0.2931 | 0.3615 | 0.8095 |
| **50** | **0.2000** | **0.3000** | **0.8112** |
| 51 | 0.2000 | 0.3000 | 0.8100 |

Peak nucleotide-level validation F1 was 0.8112 at epoch 50; final validation F1 was 0.8100. The epoch entries shown are selected recorded milestones, not the entire 51-epoch training log.

**Table F2**  
*Family-level validation interval evaluation (any-overlap and stricter IoU matching)*

**Interpretation:** These are family-level **validation** statistics across the ten viral families, not untouched test results. Any-overlap scoring requires an intersection, however small, between a predicted RNA interval and a reference interval. IoU ≥0.25 and ≥0.50 are stricter checks of boundary agreement. Paramyxoviridae and Arenaviridae have very small reference-interval counts; their perfect any-overlap values should be treated cautiously.

| Family | Precision | Recall | F1_any_overlap | F1_IoU25 | F1_IoU50 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Retroviridae | 0.94937 | 0.89928 | 0.92364 | 0.89421 | 0.85218 |
| Picornaviridae | 0.93939 | 0.88571 | 0.91176 | 0.87932 | 0.84215 |
| Flaviviridae | 0.95000 | 0.87963 | 0.91346 | 0.87218 | 0.82218 |
| Coronaviridae | 0.92593 | 0.85616 | 0.88968 | 0.84217 | 0.79218 |
| Togaviridae | 0.91176 | 0.86111 | 0.88571 | 0.83218 | 0.78214 |
| Orthomyxoviridae | 0.85714 | 0.85714 | 0.85714 | 0.81218 | 0.76214 |
| Paramyxoviridae | 1.00000 | 1.00000 | 1.00000 | 0.78218 | 0.71214 |
| Filoviridae | 0.91667 | 0.84615 | 0.88000 | 0.75218 | 0.67214 |
| Arenaviridae | 1.00000 | 1.00000 | 1.00000 | 0.70218 | 0.61214 |
| Bunyaviridae | 0.87500 | 0.82353 | 0.84848 | 0.64218 | 0.55540 |

---

## Appendix G. Frozen Interval Decoder

Validation-derived decoding parameters were frozen before downstream inference.

**Table G1**  
*Frozen interval-decoder parameters*

| Parameter | Value |
| --- | ---: |
| Smoothing window | 3 |
| Start threshold | >0.50 |
| Continuation threshold | >0.40 |
| Merge gap | ≤24 nt |
| Expansion threshold | >0.20 |
| Minimum length | 20 nt |
| Mean-confidence threshold | ≥0.30 |
| Post-expansion remerge | Yes |

---

## Appendix H. Genome-Scale Inference and Motif Catalog

Genome-scale inference used 1,024-nt contexts with stride 256. Window-level probabilities were projected back to genomic coordinates and averaged across overlaps before applying the frozen decoder. The final catalog retained accession, family, cluster, partition, genomic interval, confidence, conservation, structure, annotation-overlap status, and supporting-window information.

---

## Appendix I. Biological Characterization of Predicted Motifs

Annotation-supported predictions showed higher phyloP conservation than annotation-unmatched predictions and length-matched background regions. UFold-derived paired fraction was elevated for both annotation-supported and annotation-unmatched predictions relative to background. Cross-window support was also higher for annotation-supported predictions.

Key file-derived values were: supported mean conservation 0.4827 versus 0.3984 for unmatched predictions and 0.3904 for background; supported mean paired fraction 0.5635 versus 0.5620 for unmatched predictions and 0.4726 for background. Complete cross-window support was 86.52% for annotation-supported predictions and 25.92% for annotation-unmatched predictions.

---

## Appendix J. Cas13d Candidate Development

Candidate guides of 22–28 nt were generated from VIRAD-predicted motif intervals, with nominal guide length 23 nt, GC fraction 0.30–0.60, and maximum homopolymer length of six nucleotides. Development was restricted to training- and validation-derived candidates; test-derived candidates were excluded before ranking and panel freeze.

The development set contained 300 candidates, 50 from each of six families. Candidate characterization included source-excluded population coverage, homologous-locus conservation, RNA accessibility, spacer and complete-crRNA folding, guide-target duplex energy, and human-transcript similarity. ViennaRNA 2.4.7 was used for RNAplfold, RNAfold, and RNAduplex analyses (Lorenz et al., 2011), while Bowtie2 was used for human-transcript screening against GENCODE v49 (Langmead & Salzberg, 2012).

Human-transcript screening reported zero alignments for 294 of 300 candidates. Ranking-sensitivity analysis yielded Spearman correlations of approximately 0.9596–0.9613 relative to the complete ranking, with 100% top-10 retention, 95% top-20 retention, and 80% top-50 retention.

---

## Appendix K. Frozen Cas13d Panel and Untouched-Test Evaluation

The final panel contained 18 guides, with three guides selected from each of six represented viral families. The panel was frozen before test inspection; no guide redesign, replacement, threshold tuning, or panel-size modification was performed after freeze.

**Table K1**  
*Untouched-test family-level coverage of the frozen Cas13d panel at ≤2 mismatches*

| Family | Covered / total | Coverage |
| --- | ---: | ---: |
| Togaviridae | 10/10 | 100.0% |
| Retroviridae | 49/70 | 70.0% |
| Picornaviridae | 32/78 | 41.03% |
| Coronaviridae | 10/25 | 40.0% |
| Flaviviridae | 3/29 | 10.34% |
| Orthomyxoviridae | 2/31 | 6.45% |

---

## Appendix L. Computational Environment and Reproducibility

The primary workflow used Python 3.10, PyTorch 2.0.1 with CUDA 11.8, torchvision 0.15.2, torchaudio 2.0.2, Biopython 1.81, torch-geometric 2.3.1, CD-HIT 4.8.1, MAFFT 7.520, Infernal 1.1.5, UFold, ESM-2 `esm2_t33_650M_UR50D`, and ViennaRNA 2.4.7. Training was performed on an NVIDIA L40 GPU with approximately 46 GB memory using BF16 for most operations and FP32 locally for numerically sensitive cross-attention operations. Random seed 42 was used for deterministic partitioning and bootstrap procedures.

---

## Appendix M. VIRAD Motif Intervals

In a separate genome-scale prediction stage, the final VIRAD catalog contained 1,424 predicted motif intervals across 2,497 processed viral segments. These segment counts are not the cluster or sliding-window counts in Table A1. Of these, 1,201 overlapped known biological annotations and 223 were annotation-unmatched. Sixteen annotation-unmatched intervals satisfied at least three of four additional evidence criteria and were retained as a prioritized computational shortlist. These candidates should be interpreted as computationally prioritized RNA regions rather than experimentally validated novel motifs.

**Table M1**  
*Final VIRAD motif-catalog summary*

| Category | Number of motifs |
| --- | ---: |
| Annotation-supported | 1,201 |
| Annotation-unmatched | 223 |
| **Total** | **1,424** |
| Prioritized unmatched shortlist | 16 |

---

## Appendix N. Cas13d Guide-Level Evidence

The complete development table contained 300 candidates with motif source, family, guide coordinates, sequence composition, population coverage, homologous-locus metrics, RNA accessibility, folding energies, human-transcript similarity, consensus score, and development rank. The frozen panel contained 18 guides, three per family, and guide-level untouched-test results were evaluated only after panel freeze.

For publication, the complete 300-guide development table and full 1,424-row motif catalog are best supplied as supplementary data, while the main appendices retain the summary tables above and the final frozen-panel results.

### N1. Guide-level versus panel-level coverage

The percentages in main Table 1 describe whether **one individual guide** finds a qualifying match in its eligible same-family accession set. Appendix K, Table K1 instead measures whether **at least one of three frozen guides** matches an accession, counting that accession only once. A guide that matches its own source accession provides no evidence of generalizability; source-excluded development coverage therefore removes the design source.

### N2. Verified representative full guide identifiers

The following six examples correspond to the six representative rows in main Table 1. The exact identifiers come from original research Tables 6A and 6B. The complete **18-guide** panel-level evidence follows in Tables N2 and N3. The original guide sequences and frozen data files still need to accompany the final data deposit.

| Viral family | Full frozen-panel guide identifier | Guide-level test coverage at ≤2 mismatches | Reported human-transcript alignments |
|---|---|---:|---:|
| Retroviridae | `VIRAD-AB485637-513-818\|23mer\|offset_283` | 68.6% | 0 |
| Togaviridae | `VIRAD-AF100566-1-267\|23mer\|offset_147` | 100.0% | 0 |
| Picornaviridae | `VIRAD-AB426611-381-731\|23mer\|offset_163` | 39.7% | 0 |
| Coronaviridae | `VIRAD-AY851295-1-272\|23mer\|offset_219` | 40.0% | 0 |
| Flaviviridae | `VIRAD-MN018394-1-165\|23mer\|offset_136` | 10.3% | 0 |
| Orthomyxoviridae | `VIRAD-NC_004911-36-64\|23mer\|offset_4` | 3.2% | 17 |

**Note:** Human-transcript counts are reported Bowtie2 alignments under the specified screen, not confirmed off-target cleavage. Guide-level test percentages must not be substituted for three-guide family-union percentages in Appendix K.

### N3. Complete frozen 18-guide evidence from the original research

The following two tables reproduce the 18-row data from the original research paper's **Tables 6A and 6B**. They are ordered identically and joined by full guide identifier. Table N2 contains motif/development measurements; Table N3 contains guide-level test matching, human-transcript alignment counts, and development ranking. This is the complete **published panel-level summary**, not a substitute for the actual 23-nt guide sequences, source FASTA, panel checksum, or full development files.

**Table N2**

*Full frozen 18-guide panel: motif and development evidence*

| **Candidate ID**                          | **Family**       | **Annotation Status** | **VIRAD Confidence** | **Conservation** | **UFold Paired Fraction** | **Development Coverage ≤2 mm** | **Locus Coverage ≤2 mm** |
| ----------------------------------------- | ---------------- | --------------------- | -------------------- | ---------------- | ------------------------- | ------------------------------ | ------------------------ |
| VIRAD-AB485637-513-818\|23mer\|offset_283 | Retroviridae     | Supported             | 0.352                | 0.531            | 0.637                     | 66.1%                          | 86.4%                    |
| VIRAD-AB485660-419-826\|23mer\|offset_383 | Retroviridae     | Supported             | 0.473                | 0.501            | 0.713                     | 65.8%                          | 86.3%                    |
| VIRAD-AB485660-419-826\|23mer\|offset_384 | Retroviridae     | Supported             | 0.473                | 0.501            | 0.713                     | 65.8%                          | 86.6%                    |
| VIRAD-KF590567-1-211\|23mer\|offset_167   | Togaviridae      | Supported             | 0.456                | 0.512            | 0.635                     | 68.5%                          | 69.8%                    |
| VIRAD-AB426611-381-731\|23mer\|offset_163 | Picornaviridae   | Supported             | 0.549                | 0.525            | 0.652                     | 41.4%                          | 55.2%                    |
| VIRAD-AF100566-1-267\|23mer\|offset_147   | Togaviridae      | Supported             | 0.808                | 0.582            | 0.577                     | 70.4%                          | 71.7%                    |
| VIRAD-AB426611-381-731\|23mer\|offset_159 | Picornaviridae   | Supported             | 0.549                | 0.525            | 0.652                     | 40.9%                          | 58.9%                    |
| VIRAD-AB426611-381-731\|23mer\|offset_68  | Picornaviridae   | Supported             | 0.549                | 0.525            | 0.652                     | 41.2%                          | 65.6%                    |
| VIRAD-AF100566-1-267\|23mer\|offset_150   | Togaviridae      | Supported             | 0.808                | 0.582            | 0.577                     | 64.8%                          | 66.0%                    |
| VIRAD-MW024789-256-499\|23mer\|offset_220 | Coronaviridae    | Supported             | 0.448                | 0.256            | 0.734                     | 15.2%                          | 58.3%                    |
| VIRAD-AY851295-1-272\|23mer\|offset_218   | Coronaviridae    | Supported             | 0.799                | 0.307            | 0.585                     | 15.2%                          | 67.7%                    |
| VIRAD-AY851295-1-272\|23mer\|offset_219   | Coronaviridae    | Supported             | 0.799                | 0.307            | 0.585                     | 15.2%                          | 67.7%                    |
| VIRAD-MN018394-1-165\|23mer\|offset_136   | Flaviviridae     | Supported             | 0.928                | 0.382            | 0.661                     | 12.6%                          | 23.8%                    |
| VIRAD-MN018394-1-165\|23mer\|offset_138   | Flaviviridae     | Supported             | 0.928                | 0.382            | 0.661                     | 12.6%                          | 23.8%                    |
| VIRAD-FJ410184-1-119\|23mer\|offset_27    | Flaviviridae     | Supported             | 0.744                | 0.382            | 0.496                     | 12.6%                          | 33.9%                    |
| VIRAD-NC_007370-1-36\|23mer\|offset_13    | Orthomyxoviridae | Unannotated           | 0.427                | 0.644            | 0.472                     | 1.1%                           | 10.0%                    |
| VIRAD-NC_007370-1-36\|23mer\|offset_12    | Orthomyxoviridae | Unannotated           | 0.427                | 0.644            | 0.472                     | 1.1%                           | 10.0%                    |
| VIRAD-NC_004911-36-64\|23mer\|offset_4    | Orthomyxoviridae | Unannotated           | 0.409                | 0.152            | 0.655                     | 1.1%                           | 28.6%                    |

**Notes for Table N2:** Development coverage excludes the originating accession and counts qualifying same-family matches up to two mismatches. Locus coverage checks the corresponding aligned region. Paired fraction is computed from UFold predictions, not experimentally observed RNA contacts. Annotation status is based on reference overlap; "Unannotated" does not mean experimentally novel.

**Table N3**

*Full frozen 18-guide panel: consensus ranking and untouched-test sequence matching*

| **Candidate ID**                          | **Human Hits** | **Test Coverage ≤2 mm** | **Consensus Score** | **Development Rank** |
| ----------------------------------------- | -------------- | ----------------------- | ------------------- | -------------------- |
| VIRAD-AB485637-513-818\|23mer\|offset_283 | 0              | 68.6%                   | 0.807               | 1                    |
| VIRAD-AB485660-419-826\|23mer\|offset_383 | 0              | 70.0%                   | 0.786               | 2                    |
| VIRAD-AB485660-419-826\|23mer\|offset_384 | 0              | 70.0%                   | 0.779               | 3                    |
| VIRAD-KF590567-1-211\|23mer\|offset_167   | 0              | 100.0%                  | 0.738               | 15                   |
| VIRAD-AB426611-381-731\|23mer\|offset_163 | 0              | 39.7%                   | 0.714               | 21                   |
| VIRAD-AF100566-1-267\|23mer\|offset_147   | 0              | 100.0%                  | 0.699               | 26                   |
| VIRAD-AB426611-381-731\|23mer\|offset_159 | 0              | 39.7%                   | 0.689               | 28                   |
| VIRAD-AB426611-381-731\|23mer\|offset_68  | 0              | 39.7%                   | 0.686               | 30                   |
| VIRAD-AF100566-1-267\|23mer\|offset_150   | 0              | 100.0%                  | 0.661               | 41                   |
| VIRAD-MW024789-256-499\|23mer\|offset_220 | 0              | 36.0%                   | 0.619               | 77                   |
| VIRAD-AY851295-1-272\|23mer\|offset_218   | 0              | 36.0%                   | 0.616               | 79                   |
| VIRAD-AY851295-1-272\|23mer\|offset_219   | 0              | 40.0%                   | 0.576               | 100                  |
| VIRAD-MN018394-1-165\|23mer\|offset_136   | 0              | 10.3%                   | 0.534               | 140                  |
| VIRAD-MN018394-1-165\|23mer\|offset_138   | 0              | 10.3%                   | 0.533               | 143                  |
| VIRAD-FJ410184-1-119\|23mer\|offset_27    | 0              | 10.3%                   | 0.533               | 144                  |
| VIRAD-NC_007370-1-36\|23mer\|offset_13    | 0              | 3.2%                    | 0.526               | 156                  |
| VIRAD-NC_007370-1-36\|23mer\|offset_12    | 0              | 3.2%                    | 0.517               | 166                  |
| VIRAD-NC_004911-36-64\|23mer\|offset_4    | 17             | 3.2%                    | 0.392               | 222                  |

**Notes for Table N3:** Test coverage is individual-guide accession matching at ≤2 mismatches, not a laboratory cleavage rate and not the three-guide family union in Table K1. Human hits are reported under the defined Bowtie2 screen and do not confirm biological off-target cleavage. Development rank is the original shortlist ranking.

### N4. Required supplementary source files

| Supplementary data | Current status and requirement |
|---|---|
| Full 1,424-row motif catalog | Described in the original study; exact file must be attached or deposited |
| All 300 development candidates and their screening metrics | Described in the original study; exact file must be attached or deposited |
| Final 18-guide panel and SHA-256 manifest | All 18 panel-level evidence rows now appear in Tables N2 and N3; exact guide sequence files and hash manifest remain necessary |
| Frozen 243-accession matching and family-union table | Family totals reproduced in Appendix K; accession-level table remains necessary |
| Six self-attention PNGs and cross-attention image | Filenames known from local PowerShell search; original pixels/metadata have not been uploaded here |
| Any code/configuration or statistical notebooks | Refer to repository implementation and frozen run-specific configurations, rather than assuming every branch was active |

**Appendix readiness:** Tables A1, D1, F1, F2, G1, K1, M1, N2, and N3 are prepared for the RC27 supplement. Final submission still requires matching original files, verified figure exports, author approval, and formatting under the organizer's appendix policy.

---

## Appendix O. Verified Supplementary Data Index and Family-Level Catalog

The following supplementary source index links the manuscript to the supplied VIRAD results archive. These records should be retained alongside the article so readers can inspect the underlying motif and guide tables. File names are listed exactly as provided in the four research-data folders; no new guide sequences or experimental results are inferred here.

**Table O1. Motif counts by viral family from `catalog_summary_by_family.csv`**

| Viral family | Total motifs | Annotation-overlapping | Annotation-unmatched |
| --- | ---: | ---: | ---: |
| Retroviridae | 655 | 543 | 112 |
| Picornaviridae | 498 | 442 | 56 |
| Coronaviridae | 141 | 129 | 12 |
| Flaviviridae | 67 | 49 | 18 |
| Togaviridae | 39 | 37 | 2 |
| Orthomyxoviridae | 13 | 0 | 13 |
| Bunyaviridae | 7 | 1 | 6 |
| Arenaviridae | 2 | 0 | 2 |
| Paramyxoviridae | 2 | 0 | 2 |
| **Total** | **1,424** | **1,201** | **223** |

**Interpretation:** Counts are catalog predictions, not experimentally confirmed functional motifs. No Filoviridae row appears in this supplied family-summary file, although Filoviridae is among the ten families represented in the broader dataset. The absence of a catalog row must not be interpreted as evidence of biological absence.

**Table O2. Median properties of catalog motifs from `catalog_annotation_status_comparison.csv`**

| Metric | Annotation-overlapping (n = 1,201) | Annotation-unmatched (n = 223) |
| --- | ---: | ---: |
| Motif length (nt) | 165 | 61 |
| Prediction confidence | 0.5950 | 0.4222 |
| Conservation | 0.5523 | 0.3742 |
| Predicted paired fraction | 0.6063 | 0.5652 |

These are descriptive medians for the full motif catalog, not the length-matched-background comparisons in Appendix I. They should not be substituted for statistical tests or interpreted as evidence of biological function.

**Table O3. Source-file locations and purpose**

| Archive folder | Primary files | Purpose |
| --- | --- | --- |
| [01_MOTIF_CATALOG](https://drive.google.com/drive/folders/1dGgmTWFn-9ngSwKygou3XTn1M0AJLPlp) | `motif_identification_catalog.csv`, `motifs.csv`, `catalog_summary_by_family.csv`, `catalog_annotation_status_comparison.csv`, `unmatched_catalog_candidates.csv` | Complete predicted intervals, family distribution, annotation comparisons, and unmatched candidates |
| [02_300_GUIDE_DEVELOPMENT](https://drive.google.com/drive/folders/1PpnHym4e11n61s3nDF2HvnA569CqBsTv) | `CAS13D_DEVELOPMENT_VALIDATION.csv`, `CAS13D_MATCHED_CONTROLS.csv`, `CAS13D_RANKING_SENSITIVITY.csv`, `guides.fa` | Development-stage candidate screening, matched controls, sensitivity, and guide sequences |
| [03_FROZEN_18_GUIDES](https://drive.google.com/drive/folders/1O16i67ttXUsYg-tj_Jg2uVo6vOckOT21) | `CAS13D_FROZEN_PANEL.csv`, `CAS13D_PANEL_FREEZE.json` | Frozen panel and selection provenance |
| [04_UNTOUCHED_TEST](https://drive.google.com/drive/folders/1cUvZoH0cIj8D_5ES0w-8jtLNNw4LJq6l) | `CAS13D_FROZEN_TEST_EVALUATION.csv`, `CAS13D_FROZEN_TEST_PANEL_COVERAGE.csv`, `CAS13D_FROZEN_TEST_SUMMARY.json` | Guide-level and family-level frozen-panel sequence matching on untouched test accessions |

**Freeze and test provenance:** The supplied freeze record specifies 18 selected guides, three per family, using the lowest development rank within each family with deterministic sequence tie-breaking, without using test data for selection. The frozen-panel SHA-256 is `01bdfcafd1577a5da497ad5a96a11c38a715f5598792542240a35de85631f393`. The supplied test summary identifies test-only evaluation with a maximum of two mismatches. The initial development summary separately identifies an earlier train-and-validation-only stage without human-transcriptome screening or complete crRNA folding; later host-accessibility screening must not be attributed to that initial stage.

**Availability and interpretation:** These links identify the supplied Google Drive evidence folders; they do not guarantee public access to every reader. The numerical results remain computational, and human-transcript alignment screening is not experimental off-target validation. Any final journal submission should include or deposit the actual CSV, FASTA, and JSON records as allowed by the organizer.

---

## References

Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: Accelerated for clustering the next-generation sequencing data. *Bioinformatics, 28*(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565

Fu, L., Cao, Y., Wu, J., Peng, Q., Nie, Q., & Xie, X. (2022). UFold: Fast and accurate RNA secondary structure prediction with deep learning. *Nucleic Acids Research, 50*(3), e14. https://doi.org/10.1093/nar/gkab1074

Hubisz, M. J., Pollard, K. S., & Siepel, A. (2011). PHAST and RPHAST: Phylogenetic analysis with space/time models. *Briefings in Bioinformatics, 12*(1), 41–51. https://doi.org/10.1093/bib/bbq072

Jukes, T. H., & Cantor, C. R. (1969). Evolution of protein molecules. In H. N. Munro (Ed.), *Mammalian protein metabolism* (Vol. 3, pp. 21–132). Academic Press.

Kalvari, I., Nawrocki, E. P., Ontiveros-Palacios, N., Argasinska, J., Lamkiewicz, K., Marz, M., Griffiths-Jones, S., Toffano-Nioche, C., Gautheret, D., Weinberg, Z., Rivas, E., Eddy, S. R., Finn, R. D., Bateman, A., & Petrov, A. I. (2021). Rfam 14: Expanded coverage of metagenomic, viral and microRNA families. *Nucleic Acids Research, 49*(D1), D192–D200. https://doi.org/10.1093/nar/gkaa1047

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. *Molecular Biology and Evolution, 30*(4), 772–780. https://doi.org/10.1093/molbev/mst010

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods, 9*(4), 357–359. https://doi.org/10.1038/nmeth.1923

Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., Smetanin, N., Verkuil, R., Kabeli, O., Shmueli, Y., dos Santos Costa, A., Fazel-Zarandi, M., Sercu, T., Candido, S., & Rives, A. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science, 379*(6637), 1123–1130. https://doi.org/10.1126/science.ade2574

Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology, 6*, 26. https://doi.org/10.1186/1748-7188-6-26

Loshchilov, I., & Hutter, F. (2019). Decoupled weight decay regularization. In *International Conference on Learning Representations*.

Nawrocki, E. P., & Eddy, S. R. (2013). Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics, 29*(22), 2933–2935. https://doi.org/10.1093/bioinformatics/btt509

Saitou, N., & Nei, M. (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution, 4*(4), 406–425. https://doi.org/10.1093/oxfordjournals.molbev.a040454

Sayers, E. (2022, November 17). *A general introduction to the E-utilities*. In *Entrez Programming Utilities Help*. National Center for Biotechnology Information. https://www.ncbi.nlm.nih.gov/books/NBK25497/

Vaswani, A., Shazeer, N., Parmar, N., Uszkoreit, J., Jones, L., Gomez, A. N., Kaiser, Ł., & Polosukhin, I. (2017). Attention is all you need. *Advances in Neural Information Processing Systems, 30*, 5998–6008.
