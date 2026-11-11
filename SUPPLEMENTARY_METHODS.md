# VIRAD Supplementary Methods and Reproducibility Notes


## S1. Scope and experimental boundary

VIRAD is a **computational** study of RNA motif prediction and downstream Cas13d guide prioritization. The ten-family motif-prediction scope is distinct from the **six-family** guide-development and test scope. A predicted motif is not an experimentally confirmed functional regulatory element. Model-confidence scores, UFold structure contacts, ViennaRNA energetics, sequence-match coverage, and Bowtie2 alignments must not be presented as measured Cas13d cleavage, safety, or in-vivo antiviral activity.

## S2. Sequence acquisition, annotation and coordinate normalization

1. Query NCBI Entrez nuccore for each viral family, prioritizing RefSeq and including constrained additional complete-genome records. Preserve raw FASTA/GFF3, accession/version, family, source, sequence length, and provenance (Sayers, 2022).
2. Normalize accession joins without silently merging different biological records. Preserve the original identifiers in all derived tables.
3. Combine NCBI RNA-relevant GFF3 intervals with Rfam/Infernal covariance-model evidence. In `cmscan`, the study used `--cut_ga --rfam --nohmmer` and a structured tabular output (Kalvari et al., 2021; Nawrocki & Eddy, 2013).
4. Label motif presence per nucleotide as binary, and use the 87-class motif-type target only when a valid specific label exists. Missing motif-type labels use ignore index −1; incompatible overlapping labels are not arbitrarily resolved.
5. Convert external one-based inclusive genomic coordinates into **zero-based half-open** intervals `[start, end)`, with interval length `end - start`. Distinguish alignment columns from ungapped genomic positions.

**Failure rule:** A mismatch among accession, window, annotation, or modality coordinate systems invalidates the example rather than triggering an undocumented realignment.

## S3. Sequence-quality thresholds and cluster-disjoint partitioning

| Stage | Frozen criterion |
|---|---|
| Whole genome | Length at least 500 nt |
| Whole genome ambiguity | N fraction no more than 0.10 |
| Local 512-nt window | N fraction no more than 0.20 |
| Local 512-nt window | Gap fraction no more than 0.20 |
| CD-HIT-EST clustering | 95% identity; 90% coverage of both aligned sequences |
| Stratification | Viral family and motif-label presence |
| Cluster assignment | About 70% train, 15% validation, 15% test, using seed 42 |
| Disjointness | No cluster may occur in more than one partition |
| Window construction | 512 nt, 256-nt stride; no primary trailing-window padding |

Within each stratum containing `n` clusters, training allocation was `round(0.70*n)`, validation `round(0.15*n)`, and the remainder test, using Python's standard rounding semantics. Individual sequences and all derived windows inherit their source cluster's assignment.

**Frozen observed data:** train 59,630 windows/1,724 clusters; validation 12,852 windows/372 clusters; test 13,209 windows/375 clusters. Totals: **85,691 windows and 2,471 clusters**, with zero shared cluster identifiers.

## S4. Conservation model fitted only to training records

Within each family, MAFFT aligns training genomes (Katoh & Standley, 2013). Pairwise comparisons count mismatches only where *both* aligned rows contain canonical A/C/G/T. For `N_ab` comparable sites, the mismatch fraction is `p_ab = mismatches/N_ab`. A Jukes–Cantor transformed distance is `d = −(3/4) ln(1 − 4p/3)` with the numerical cap `p* = min(p, 0.75 − 10^−6)` to keep the logarithm defined. Sequence pairs with no canonical overlap are unsupported, not assigned an invented distance. Deterministic neighbor joining produces a family tree. PHAST `phyloFit` uses a REV model, and `phyloP` uses likelihood-ratio `CONACC` with reference index zero (Hubisz et al., 2011).

Validation and test sequences are projected into the **training-defined** evolutionary background without refitting the training alignment, tree, or substitution model. The mapped conservation feature is bounded as `min(1, max(0, x))`. This bounded input is **not** a probability of motif function.

## S5. Per-window modalities and data-validation contract

| Component | Input shape for length L | Interpretation |
|---|---|---|
| Nucleotide sequence | `[B,5,L]` | One-hot A/U/G/C/N |
| PhyloP conservation | `[B,L]` | Training-derived evolutionary score |
| UFold structure | `[B,L,L]` | Symmetric binary contact map; **not calibrated probabilities** |
| ESM-2 host context | `[B,1280]` | Global associated host-protein representation |
| Attention mask | `[B,L]` | Valid nucleotide positions |
| Motif presence | `[B,L]` | Nucleotide-level binary target |
| Motif type | `[B,L]` | Class 0–86 or ignore index −1 |
| Shared model representation | `[B,L,768]` | Fused contextual state |

The UFold contact map is accepted only if its cache key, dimensionality, nucleotide-origin mapping, finiteness, symmetry and binary values are consistent. An all-zero valid map is not an automatic error. The model input's host representation is built with ESM-2 `esm2_t33_650M_UR50D`; it is projected from 1,280 to 768 dimensions and cannot be interpreted as observed interaction at each viral nucleotide (Fu et al., 2022; Lin et al., 2023).

## S6. Network design and curriculum optimization

The main network comprises **68,819,416 parameters**, six Viroformer blocks, hidden width 768, 12 self-attention heads, and feed-forward width 3,072. Biological signals enter learned modality-fusion gates, conservation and structure attention biases, a six-head sequence-to-structure cross-attention path, host cross-attention, and one symmetry-enforced structural pair-refinement pass. Dropout and stochastic depth were each 0.10. Structure-related attention weights are learned tendencies rather than a deterministic claim of RNA pairing.

| Stage | Maximum epochs | Enabled information | Initial learning rate |
|---|---:|---|---:|
| 1 | 2 | Sequence | 1 × 10⁻⁴ |
| 2 | 5 | Sequence, conservation, RNA structure | 1 × 10⁻⁴ |
| 3 | 5 | Above plus host context | 5 × 10⁻⁵ |
| 4 | 43 | Full multimodal/multitask optimization | 1 × 10⁻⁵ |

Configured maximum: **55 epochs**; actual final continuation: **51 epochs and 95,064 optimizer steps**. AdamW uses β₁=0.9, β₂=0.999, weight decay 0.01; global gradient clipping uses 1.0. A micro-batch of 4 with eight accumulated batches gives effective batch 32. The exponential-moving-average teacher has decay 0.999. Warmup and cosine learning-rate decay approach 10⁻⁶. On the NVIDIA L40 system, BF16 mixed precision was combined with local FP32 execution for numerically sensitive cross-attention operations.

**Model-selection rule:** The selected checkpoint changes on strict validation-loss improvement. **Early-stopping rule:** Stage-four validation F1, patience 8, minimum global epoch 20. These are distinct mechanisms. No held-out test examples contribute to gradients or selection decisions.

## S7. Validation-controlled interval decoding and whole-genome mapping

The nucleotide probability track is converted to genomic intervals by a deterministic decoder selected on validation data:

1. Smooth with a width-3 same-mode convolution.
2. Create seeds where smoothed probability **>0.50**.
3. Continue through connected positions with probability **>0.40**.
4. Merge neighboring regions separated by **≤24 nt**.
5. Expand boundaries while adjacent smoothed probability **>0.20**.
6. Exclude intervals shorter than **20 nt**.
7. Exclude intervals whose mean smoothed probability is **<0.30**.
8. Remerge touching or overlapping expanded intervals.

Threshold inequalities and the **order of steps** must remain unchanged. Whole-genome inference uses overlapping 1,024-nt windows at 256-nt stride, projects local scores to source coordinates, aggregates overlapping evidence, and applies the frozen decoder. The catalog retains confidence, provenance, interval coordinates, annotations, and structure/conservation summaries.

## S8. Motif biological characterization

* Compare each predicted interval against RNA-relevant reference annotation, with the two categories **annotation-supported** and **annotation-unmatched**. Unmatched does not mean novel or functional.
* Compare nucleotide conservation among existing annotated motifs, supported predictions, unmatched predictions, and genomic background.
* Examine UFold-predicted pairedness and window-overlap reproducibility, ensuring background comparisons use consistent biological units.
* Preserve non-significant results: annotation-unmatched predictions did **not** show statistically significant overall conservation enrichment against background in the reported analyses, although they exhibited greater predicted RNA pairing.

**Observed final catalog:** 1,424 candidate motifs from 2,497 viral segments; 1,201 (84.3%) annotation-supported; 223 (15.7%) annotation-unmatched.

## S9. Development guide enumeration and data firewall

Eligible targets come **only from training/validation-derived motif intervals**. Enumerate **22–28-nt** targets with nominal **23-nt** design length, GC fraction 0.30–0.60, and maximum homopolymer length 6. Save source accession, family, motif interval ID, genomic coordinate and target. When segmented accession sequences are concatenated for searching, insert 50 N bases between segments. Balanced development guide shortlist: **300 candidates** across Coronaviridae, Flaviviridae, Orthomyxoviridae, Picornaviridae, Retroviridae, and Togaviridae.

The crRNA spacer must be complementary to its RNA target; the VIRAD interval is the *target* region, not the same as the full crRNA scaffold. Guide design, filtering, scoring, ranking, and panel freeze do not inspect test-derived candidates.

## S10. Coverage, homologous-locus comparison, and mismatch criteria

* **Source-excluded whole-accession occurrence:** Remove each guide's design-source accession before searching all other same-family development accessions. Search guide and reverse-complement orientations across equal-length candidate windows. Record minimum Hamming distance; thresholds are 0, ≤1, and ≤2 nucleotide differences. Nominal 23-mers are accelerated by an 8+8+7 seed strategy but full-length Hamming verification determines matches.
* **Locus conservation:** Align development accessions within families; map guide coordinates to homologous alignment positions. Exclude guide source, uninformative and gapped loci. Summarize per-position mismatch rates, haplotype diversity, major haplotypes, and Shannon nucleotide entropy `H_i = −Σ_b p_i,b log2(p_i,b)`. Whole-genome occurrence and homologous-locus support answer different questions.
* **Development family-alignment populations:** Coronaviridae 139; Flaviviridae 152; Orthomyxoviridae 175; Picornaviridae 431; Retroviridae 396; Togaviridae 55 accessions, as reported in the study documentation.

## S11. Structural and host-transcript guide checks

| Evidence | Computational method and meaning |
|---|---|
| Target accessibility | RNAplfold local probability that the *entire guide-length interval* is unpaired |
| Guide self-folding | RNAfold minimum free energy for guide-complementary spacer |
| Complete-crRNA folding | RNAfold on direct repeat + spacer |
| Guide–target duplex | RNAduplex interaction energy between antisense guide and target |
| Human-transcript similarity | Bowtie2 alignments against GENCODE v49 transcriptome |
| Matched local controls | Same accession, within ±500 nt, equal length, GC difference ≤0.05, no guide overlap |

ViennaRNA settings: `-W min(150,n)`, `-L min(100,max(1,n−1))`, `-u 28`. RNAplfold `*_lunp` values were read for the complete candidate interval using offset and interval-length indexing defined by the original pipeline. The full-crRNA direct repeat used in the computational model was `AACCCCTACCAACTGGTCGGGGTTTGAAAC`; its sequence alone does not validate Cas13d loading or imply ortholog-independent functionality.

Bowtie2 used end-to-end alignment, `-N 1`, `-L 15`, and `-k 100`. Thus a guide's count of reported hits is capped, and no-hit results do **not** guarantee off-target safety. The observed screen reported no transcript alignments for **294/300 candidates** under that configuration.

## S12. Guide rank integration and sensitivity

Each available guide metric was converted to a direction-aware percentile rank. Metrics with larger preferred values included exact and one-mismatch development coverage, exact homologous-locus support, accessibility, motif conservation, and motif prediction confidence. Metrics with smaller preferred values included detected human-transcript alignments and locus mean entropy. Tied values received mean percentile rank; missing optional evidence was omitted rather than imputed. The consensus score equals the mean of a guide's available transformed percentiles. Tie-breaking used descending exact coverage, descending one-mismatch coverage, and then target-sequence order. Scores are **relative ranking indices**, not calibrated Cas13d activity probabilities.

Sensitivity analyses recalculated ranks under predetermined available evidence subsets. The reported ten highest-ranked guides remained among the top ten across scoring methods tested, but this does not establish universal robustness to all possible weighting choices.

## S13. Frozen-panel integrity and untouched test

Select the three highest-ranked development guides from **each of the six** guide families, using prespecified deterministic ordering. Freeze panel size, selection rules and 18 guide sequences **before** test access. Compute SHA-256 hashes of the development table and frozen panel. Verify the frozen panel's observed hash against the freeze record; abort if changed. Only then measure coverage on **243 never-used test accessions**, separated by family:

| Family | Untouched test accessions |
|---|---:|
| Coronaviridae | 25 |
| Flaviviridae | 29 |
| Orthomyxoviridae | 31 |
| Picornaviridae | 78 |
| Retroviridae | 70 |
| Togaviridae | 10 |
| **Total** | **243** |

Test only each guide against its same-family accession set. Compute exact, ≤1- and ≤2-mismatch guide coverage and the **family union** at ≤2 mismatches: each accession contributes at most once when *any* of the family's three guides matches. No test-informed optimization, new guide selection, panel-size alteration, or threshold change is allowed.

**Final ≤2-mismatch family coverage reported:** Togaviridae 10/10 (100%); Retroviridae 49/70 (70%); Picornaviridae 41.03%; Coronaviridae 40.00%; Flaviviridae 10.34%; Orthomyxoviridae 6.45%. The different family sample sizes must appear alongside percentages. These are **in-silico sequence matches**, not validated cleavage efficiencies.

## S14. Statistical analysis

The unit of observation is chosen according to the analysis: accession, sequence, cluster, or family, not by counting overlapping windows as independent organisms. Paired and family-stratified bootstrap analyses used **10,000 draws** with seed 42; the latter resampled each family at its original size. Report 95% percentile confidence intervals from the 2.5th and 97.5th quantiles. Multiple hypothesis tests used Benjamini–Hochberg false-discovery-rate correction or Holm family-wise correction as specified for the comparison. Distinguish (a) prediction scores, (b) interval detection, (c) annotation agreement, (d) biological feature enrichment, and (e) independent guide sequence coverage; they are not interchangeable endpoints.

## S15. Documentation scope

The supplementary materials present computational methods, validation procedures, and numerical evidence through text, equations, and tables.

## S16. Supplementary documentation

This section documents the study's computational methods, model specifications, and validation boundaries. Associated numerical evidence is provided in [Supplementary Tables](APPENDICES_AND_TABLES.md), and the supporting data are available in the [VIRAD data archive](https://drive.google.com/drive/folders/1n6sDwAfD4rUDg9qRLN996hVdfkCDaJew).

## S17. Canonical Mathematical Equations (Full Set)

The equations below reproduce the labeled mathematical definitions from the detailed methodology (original numbering retained). They are supplied for scientific auditability; a few representative equations appear in the short article. Symbols are explained in the corresponding method descriptions above and in the canonical methodology. These are **computational definitions**, not experimentally measured biological laws.

**Equation 1 | Half-open interval length**

$$
L_{\mathrm{interval}} = e-s
$$

**Equation 2 | Ambiguous-base fraction**

$$
f_N=\frac{1}{L}\sum_{i=1}^{L}\mathbf{1}[x_i=N]
$$

**Equation 3a | Training-cluster allocation**

$$
n_{\mathrm{train}}=\operatorname{round}(0.70n)
$$

**Equation 3b | Validation-cluster allocation**

$$
n_{\mathrm{val}}=\operatorname{round}(0.15n)
$$

**Equation 3c | Test-cluster allocation**

$$
n_{\mathrm{test}}=n-n_{\mathrm{train}}-n_{\mathrm{val}}
$$

**Equation 4 | Cluster-disjointness constraint**

$$
C_{\mathrm{train}}\cap C_{\mathrm{val}}=C_{\mathrm{train}}\cap C_{\mathrm{test}}=C_{\mathrm{val}}\cap C_{\mathrm{test}}=\varnothing
$$

**Equation 5 | Pairwise mismatch fraction**

$$
p_{ab}=\frac{1}{N_{ab}}\sum_{k\in M_{ab}}\mathbf{1}[a_k\ne b_k],\qquad N_{ab}=|M_{ab}|
$$

**Equation 6 | Jukes-Cantor distance**

$$
d_{ab}=-\frac{3}{4}\ln\!\left(1-\frac{4p_{ab}^{*}}{3}\right)
$$

**Equation 7 | Jukes-Cantor saturation cap**

$$
p_{ab}^{*}=\min\!\left(p_{ab},\,0.75-10^{-6}\right)
$$

**Equation 8 | Conservation feature bounding**

$$
c_i=\min\!\left(1,\max(0,x_i)\right)
$$

**Equation 9 | Training-window start**

$$
s_k=256k
$$

**Equation 10 | Complete-window requirement**

$$
s_k+512\le N
$$

**Equation 11 | Binary RNA-contact representation**

$$
P\in\{0,1\}^{L\times L},\qquad P_{ij}=P_{ji}
$$

**Equation 12 | Learned modality weight**

$$
w_m=\frac{\exp(g_m)}{\sum_{r\in\mathcal{M}}\exp(g_r)}
$$

**Equation 13 | Direct multimodal fusion**

$$
E_{\mathrm{mix}}=\sum_{m\in\mathcal{M}}w_mE_m
$$

**Equation 14 | Content attention score**

$$
s_{ij}^{(h)}=\frac{\left(q_i^{(h)}\right)^{\!T}k_j^{(h)}}{\sqrt{64}}
$$

**Equation 15 | Conservation attention bias**

$$
b_{ij,\mathrm{cons}}^{(h)}=\ln\!\left(1+\max\!\left[\alpha_h(c_i+c_j),0\right]\right)
$$

**Equation 16 | Structure attention bias**

$$
b_{ij,\mathrm{struct}}^{(h)}=\beta_h P_{ij}
$$

**Equation 17 | Biologically biased attention logit**

$$
a_{ij}^{(h)}=s_{ij}^{(h)}+b_{ij,\mathrm{cons}}^{(h)}+b_{ij,\mathrm{struct}}^{(h)}
$$

**Equation 18 | Conservation auxiliary target**

$$
y_i^{\mathrm{cons}}=\mathbf{1}[c_i>0.7]
$$

**Equation 19 | Masked binary cross-entropy**

$$
\mathcal{L}_{\mathrm{BCE}}=\frac{\sum_i m_i\,\ell_{\mathrm{BCE}}(z_i,y_i)}{\max\!\left(\sum_i m_i,1\right)}
$$

**Equation 20 | Dice component**

$$
\mathcal{L}_{\mathrm{Dice}}=1-\frac{2\sum_i p_i y_i}{\max\!\left(\sum_i p_i+\sum_i y_i,10^{-8}\right)}
$$

**Equation 21 | Motif segmentation objective**

$$
\mathcal{L}_{\mathrm{motif}}=\mathcal{L}_{\mathrm{Dice}}+\mathcal{L}_{\mathrm{BCE}}
$$

**Equation 22 | Learned uncertainty weighting**

$$
\widetilde{\mathcal{L}}_k=\frac{1}{2}\exp(-2s_k)\mathcal{L}_k+s_k
$$

**Equation 23 | Grouped uncertainty-weighted objective**

$$
\mathcal{L}_{\mathrm{group}}=\sum_{k\in\mathcal{A}}\widetilde{\mathcal{L}}_k
$$

**Equation 24 | Total multitask objective**

$$
\mathcal{L}_{\mathrm{total}}=\mathcal{L}_{\mathrm{group}}+\mathcal{L}_{\mathrm{type}}+\mathcal{L}_{\mathrm{presence}}+\lambda_{\mathrm{KL}}\mathcal{L}_{\mathrm{KL}}
$$

**Equation 25 | Effective batch size**

$$
B_{\mathrm{eff}}=B_{\mathrm{micro}}N_{\mathrm{accum}}
$$

**Equation 26 | EMA teacher update**

$$
\phi_t=0.999\,\phi_{t-1}+0.001\,\theta_t
$$

**Equation 27 | Overlap-aggregated nucleotide probability**

$$
\bar{p}_i=\frac{1}{K_i}\sum_{k=1}^{K_i}p_i^{(k)},\qquad K_i>0
$$

**Equation 28 | Mean interval conservation**

$$
\bar{c}_I=\frac{1}{|I|}\sum_{i\in I}c_i
$$

**Equation 29 | Hamming distance**

$$
H(x,y)=\sum_{i=1}^{L}\mathbf{1}[x_i \ne y_i]
$$

**Equation 30 | Source-excluded accession coverage**

$$
\mathrm{Coverage}_t(g)=\frac{1}{|A_g|}\sum_{a\in A_g}\mathbf{1}\!\left[H_{\min}(g,a)\le t\right]\times100\%,
\qquad t\in\{0,1,2\}
$$

**Equation 31 | Nucleotide entropy**

$$
H_i=-\sum_{b\in\{A,C,G,T\}}p_{i,b}\log_2 p_{i,b}
$$

**Equation 32 | Available-metric percentile consensus score**

$$
S_g=\frac{1}{|M_g|}\sum_{m\in M_g}r_{m,g}
$$

**Equation 33 | Panel union-coverage set**

$$
A^{\mathrm{test}}_{\mathrm{panel},t}=\bigcup_{g\in G}A^{\mathrm{test}}_{g,t}
$$

**Equation 34 | Panel-level untouched-test coverage**

$$
\mathrm{Coverage}_{\mathrm{panel},t}=\frac{\left|A^{\mathrm{test}}_{\mathrm{panel},t}\right|}{|A_{\mathrm{test}}|}\times100\%
$$

**Equation 35 | Bootstrap paired mean difference**

$$
\Delta^{*(b)}=\frac{1}{n}\sum_{i\in S_b}(a_i-b_i)
$$

**Equation 36 | Percentile bootstrap interval**

$$
\mathrm{CI}_{0.95}=\left[Q_{0.025}(\Delta^{*}),\,Q_{0.975}(\Delta^{*})\right]
$$

## S18. Expanded chronological methods and audit checkpoints

This section spells out the operational sequence underlying Sections S2 to S14. It supplements the numerical thresholds already listed, rather than introducing different processing settings. Every step is computational. Where the repository does not establish a parameter or a sample count, it is not inferred.

### S18.1 Source retrieval, data provenance, and RNA-label construction

For each of the ten prespecified families, source nucleotide records were queried using the family name and RefSeq accession filters, followed by additional complete-genome records within family-specific retrieval limits. The process retained the untouched FASTA and annotation files before creating normalized join identifiers. Each derivative was associated with its accession, version when available, family assignment, original sequence length, and source. During the join of GFF3 intervals, Rfam covariance-model matches, evolutionary scores, structural predictions, window coordinates, and model targets, an identifier mismatch was considered a provenance problem. An accession key was not allowed to silently merge two different biological sequences.

RNA-relevant supervision was constructed independently of model inference. NCBI RNA-related intervals and Rfam/Infernal matches provided binary motif-presence supervision and motif-type labels when the annotation evidence was sufficiently specific. Incompatible motif-type annotations were not forced into one category. For unlabeled motif types, the target was masked with ignore index −1; this means the model was not penalized for guessing an unsupported motif class. GFF3 positions were converted from one-based inclusive to zero-based, half-open intervals before mapping them into local windows, keeping the original genomic coordinates in the provenance record.

### S18.2 Filtering, grouping, splitting, and window indexing

Whole sequences shorter than 500 nt and sequences exceeding 10% ambiguous bases were removed before local-window quality checks. Local candidate windows had separate limits of 20% ambiguous bases and 20% gaps. The sequence clustering stage used CD-HIT-EST at 95% identity and at least 90% aligned coverage in *each* sequence, avoiding short fragment similarities being treated as evidence of near-identical entire genomes. Cluster assignment used a deterministic random seed of 42 and stratification by viral family and presence or absence of RNA-motif labels. This assignment was completed **before sliding windows were materialized**, and every member of a cluster, together with all windows derived from it, inherited the same partition.

A 512-nt window was extracted every 256 nt for model preparation. The windows overlapped by half their length, which supplied local context for motifs crossing a single window's edge. The manifest preserved accession, cluster, partition, genomic start/end coordinate, source motif intervals, label-mask coverage, and feature provenance. The disjointness audit tested the intersection of train, validation, and test cluster-ID sets; all three pairwise intersections were empty. The frozen totals were 1,724, 372, and 375 clusters, producing 59,630, 12,852, and 13,209 windows, respectively. The fact that windows overlap within a partition does not make them statistically independent samples.

### S18.3 Evolutionary modeling without leakage

Evolutionary features were fitted separately by family using **training genomes only**. The alignment stage produced multiple-sequence alignments with MAFFT. Pairwise mismatch fractions were calculated over columns where both sequences had canonical bases, followed by Jukes–Cantor correction and deterministic neighbor joining. The family tree then provided a topology for phyloFit substitution-model fitting, and phyloP supplied conservation evidence at relevant nucleotide positions. Zero-comparable-site pairs were not assigned artificial evolutionary distances. Validation and test examples could be **scored** using the already frozen evolutionary background but could not revise its alignment, topology, or substitution parameters. This is a second firewall beyond cluster-disjoint partitioning: without it, model input features could leak information from held-out genomes.

### S18.4 Four synchronized modalities and internal representation

The nucleotide modality consisted of five one-hot channels representing A, U, G, C, and ambiguous or unknown identity. A corresponding conservation vector supplied the mapped phyloP-derived score at each nucleotide. UFold produced a binary, symmetric nucleotide-contact matrix; a value of one indicated a predicted pair in that output, **not a measured or calibrated pairing probability**. The protein-context modality used an ESM-2 embedding of a host-associated protein input, yielding one 1,280-component vector per eligible record. The detailed identity and availability of host-protein inputs must be retained in provenance records; the embedding alone does not establish direct protein-RNA binding.

The model admitted the window only if nucleotide indices, conservation values, contact-map dimensions, and masks were consistent. Expected batch shapes were [B,5,L] for sequence, [B,L] for conservation and masks, [B,L,L] for the contact matrix, and [B,1280] for host context, with L=512 for model training. Nonfinite values, mismatched accession-origin mappings, and invalid contact-map structure were not silently accepted. The projections formed 768-dimensional shared hidden states for Viroformer.

### S18.5 Six-block encoder, loss functions, and optimization

The 68,819,416-parameter Viroformer comprises six encoder blocks with 12 self-attention heads per block, width 768, and a feed-forward intermediate width of 3,072. Learned input projections and fusion gates combine nucleotide, conservation, structure, and host information. Conservation-derived biases and predicted structural contacts influence attention, and sequence-to-structure cross-attention and the host-context path provide additional representation mixing. Learned attention preferences are interpretability signals, not molecular interaction evidence.

The supervised learning target is the nucleotide-level motif-presence mask, with motif category as an available auxiliary prediction target. Training combines positionwise segmentation and overlap-sensitive losses with the auxiliary tasks allowed by the source annotations. Missing motif-type labels are masked from the corresponding auxiliary contribution rather than imputed. The four curriculum stages successively enable: (1) nucleotide sequence for at most 2 epochs; (2) conservation and predicted structure for 5 additional epochs; (3) host-associated context for 5 epochs; and (4) full multimodal and multitask training for at most 43 further epochs. The planned upper limit was 55 epochs; the logged run finished at epoch 51 and 95,064 optimizer updates.

Optimization used AdamW with β₁=0.9, β₂=0.999, weight decay 0.01, gradient norm clipping at 1.0, micro-batch size 4, and accumulation over eight micro-batches for effective batch size 32. The source-listed initial learning rates were 1 × 10⁻⁴ for stages 1 and 2, 5 × 10⁻⁵ for stage 3, and 1 × 10⁻⁵ for stage 4, with warmup and cosine decay towards 10⁻⁶. The moving-average teacher used 0.999 decay. Mixed-precision BF16 was applied when supported; numerically sensitive cross-attention operations were evaluated in FP32. The best checkpoint was updated on strict validation-loss improvement, while stage-four early stopping monitored validation F1 with patience eight after minimum global epoch twenty. This training progression is **not a controlled ablation study**, because stage duration, optimization state, and available inputs all change together.

### S18.6 From nucleotide predictions to genome-scale motif intervals

The model's immediate output is a probability for each nucleotide, not an RNA motif interval. A validation-selected decoder first smooths the probability series with a three-position convolution, seeds a motif where smoothed probability is above 0.50, continues through adjacent scores above 0.40, merges nearby fragments separated by at most 24 nt, and expands the boundaries while neighboring scores exceed 0.20. An interval must span at least 20 nt and have mean smoothed confidence of at least 0.30; overlapping or touching expanded intervals are merged again. This sequence of operations is important because a region with many high-scoring fragments can otherwise be incorrectly counted as numerous separate motifs.

Whole-genome inference differs from training-window generation: it uses 1,024-nt contexts with 256-nt stride to cover long viral records. Each window's probability output is mapped back to the unchanged genomic coordinates. For every nucleotide covered more than once, the respective window-level probabilities are averaged before the frozen decoder is applied. The resulting catalog stores coordinates, source accession, family, confidence, annotation overlap, conservation, structure evidence and cross-window support. The validation-selected thresholds are **not adjusted based on the final test accessions**.

### S18.7 Distinguishing annotation support from candidate novelty

After inference, each motif interval is intersected with reference RNA-related annotations. Overlap yields the category **annotation-supported**, while absence of overlap yields **annotation-unmatched**. This post-prediction label describes the completeness of the available annotation comparison, not biological novelty. The catalog comprised 1,424 intervals from 2,497 processed viral segments, of which 1,201 were annotation-supported and 223 annotation-unmatched.

Biological evidence is evaluated in complementary ways: conservation summarizes evolutionary constraint, the predicted paired fraction characterizes potential RNA structure, and overlapping-window support measures numerical reproducibility within local sequence context. A four-property prioritization scheme considered model confidence, evolutionary conservation, UFold-predicted paired fraction and cross-window support. **Sixteen annotation-unmatched intervals met at least three of the four source-defined criteria**. The source narrative identifies the criterion dimensions but does not give all cutoff values in this supplementary archive, so missing threshold values must be recovered from the frozen analysis rather than estimated. Neither an annotation-unmatched classification nor strong computational evidence establishes that a region has an experimentally confirmed RNA function.

### S18.8 Cas13d candidates, source exclusion, and screening populations

The guide generation stage uses **only development-derived motif intervals**, not motif predictions chosen by inspecting the untouched test accession population. For each eligible region, short target sequences of 22–28 nt were enumerated; 23 nt was the nominal length. First-pass filters required GC fraction between 0.30 and 0.60 and homopolymer length no greater than six. These project-specific numerical criteria simplify computational candidate selection and are not universal Cas13d biochemical thresholds. Each target retained its source accession, family, motif ID, interval coordinates and target sequence; the complementary crRNA spacer was treated distinctly from the direct-repeat scaffold.

The final development shortlist had 300 guides, 50 from each of six families. To estimate how common a candidate's target was outside its design source, the pipeline excluded its originating accession and searched the other development accessions from the same family. Whole-length candidate windows were checked for exact, one-mismatch and two-mismatch coverage, with the search orientations and full-length Hamming verification fixed in advance. An entire accession was counted as covered when a qualifying occurrence existed anywhere in it. This result is **not the same as conservation of the homologous genomic locus**. For locus conservation, separately aligned development-family genomes were mapped to corresponding positions, with gapped or uninformative loci excluded from supported matches and per-position entropy/haplotype summaries calculated. No development-stage homologous-locus model was refitted using the test collection.

### S18.9 Orthogonal computational guide checks and negative controls

RNAplfold estimated the probability that the full guide-target interval is unpaired under local folding assumptions. RNAfold evaluated spacer and complete-crRNA folding; RNAduplex computed guide-target hybridization energy. These quantities provide different computational proxies for target accessibility, self-structure and duplex thermodynamics. A strongly favorable duplex energy does not prove that a guide can reach a target in an infected cell. Local comparison targets were chosen from the same accession within 500 nt, required equal length and GC-fraction difference of at most 0.05, and excluded overlap with the candidate site. These matched controls help contextualize a guide's relative properties without testing antiviral function.

Human transcript screening used Bowtie2 end-to-end alignments against GENCODE v49 with seed length 15, one seed mismatch permitted and at most 100 reported hits per guide. The reported outcome that 294/300 development candidates yielded no alignment belongs strictly to those search parameters. Absence of reported matches does not mean that a guide has no biologically relevant off-target effects, and the capped hit count is not an exhaustive enumeration of all potential transcript interactions.

### S18.10 Consensus ranking, freeze integrity and independent sequence coverage

Guide metrics were standardized into percentile ranks according to whether larger or smaller values were desirable. The ranked evidence included population coverage, homologous-locus agreement, local entropy, predicted accessibility, motif conservation/confidence and human-transcript alignment burden. Missing optional measurements were omitted rather than replaced with arbitrary favorable values. The average of the available direction-aware ranks formed the consensus index, followed by deterministic tie-breaking. Sensitivity analysis re-ranked candidates under alternative available-metric sets; the reported top-ten stability applies to those tested alternatives and is not proof of robustness under all possible choices.

Exactly three top-ranked guides per family were selected **using development scores only**. The selected 18-guide panel, its ranking source and preset comparison threshold were frozen before any untouched-test assessment, and SHA-256 checks confirmed that the evaluated panel matched the freeze record. The test population comprised 243 never-used viral accessions across the same six guide families, with counts 25, 29, 31, 78, 70 and 10 for Coronaviridae, Flaviviridae, Orthomyxoviridae, Picornaviridae, Retroviridae and Togaviridae, respectively. Each guide was searched only against its same-family test accession set at exact, at most one mismatch, and at most two mismatches. Family-union coverage counted an accession only once if any of that family's three guides qualified. **No guide replacement, threshold retuning, panel-size adjustment, or family-specific selection changes were permitted after test results were inspected.** These percentages demonstrate computational target-sequence coverage, not Cas13d cleavage.

### S18.11 Statistical unit, uncertainty and failure modes

Overlapping training windows from one viral genome are correlated observations, so interval-level, accession-level, cluster-level and family-level summaries should not be interpreted as interchangeable sample sizes. Paired bootstrap comparisons and family-stratified paired resampling used 10,000 replicates with seed 42, reporting 95% percentile intervals where the source results provide them. Multiplicity correction followed the comparison-specific Benjamini–Hochberg or Holm procedures; this supplementary archive does not invent adjusted p values or interval endpoints absent from the frozen outputs. The inference controls also prevent a high any-overlap F1 from being described as precise motif-boundary recovery, attention entropy from being described as RNA contact evidence, or guide sequence matches from being described as experimentally measured antiviral efficacy.

### S18.12 Reproduction package and final source checks

A complete reproducibility package should retain the original NCBI accession manifests and FASTA/GFF3 downloads, annotation join and coordinate-normalization logs, CD-HIT-EST assignments, partition manifest and shared-cluster verification, training-only MAFFT/phyloFit/phyloP artifacts, cached UFold features and host-protein embedding provenance, model configuration and selected checkpoint, interval-decoder parameters and overlapping-window score tracks, the final motif catalog, the 300-guide development tables and matched controls, and the frozen 18-guide panel with its integrity hashes and untouched accession match table. The documented workflow and summary results inform this **recommended audit inventory**; the list does not establish that every artifact is available in the archive.

## S19. Viroformer architecture, equation-level implementation, and information flow

This section expands the study Multimodal Modeling description using the implementation modules in `src/virad/models/embeddings.py`, `gnn.py`, `fusion.py`, `attention.py`, `viroformer.py`, `structure_module.py`, `heads.py`, and `virad.py`. It is a description of the **implemented computational model**, not an independent biological validation or an assertion that every optional learning branch was active in every training run.

### S19.1 Shared coordinate system and modality projections

For a batch of (B) windows of length (L=512), the nucleotide one-hot input has shape ([B,5,L]). It is transposed into ([B,L,5]) and projected by a bias-free linear mapping into ([B,L,768]). Sinusoidal positional vectors of the same width are added before layer normalization and dropout, retaining ordered genomic position information within the window.

The conservation track has shape ([B,L]). At each position, a two-layer network maps one scalar through a 128-unit hidden layer with ReLU into width 768, followed by layer normalization. The UFold structure is a ([B,L,L]) symmetric binary adjacency matrix. By default, a three-layer graph network updates each nucleotide state using both its existing sequence-derived representation and a sum of its neighbors according to the contact matrix. At each layer the graph operation is `ReLU(W0 h + W1 (P h))`, followed by a projection and normalization at output. This supplies the structural branch with information about predicted contacts, including contacts between positions far apart along the nucleotide sequence.

The host-associated input is a single 1,280-element ESM-2-derived protein representation per example. A linear 1,280-to-768 projection and layer normalization produce a global host vector, and the vector is broadcast to every nucleotide position. Broadcasting allows host context to condition the model but does **not** create independent nucleotide-level host interaction evidence.

### S19.2 Learned four-way fusion and sequence-to-structure attention

The fusion module stores four learnable scalar gate logits, initialized equally. A softmax converts them into nonnegative weights summing to one. The four 768-dimensional modality representations are added as a gate-weighted sum. The implemented logits are shared at the fusion level, not separately computed from a different gate for each nucleotide. The source equations are:

```math
w_m = \frac{\exp(g_m)}{\sum_{r \in \mathcal M}\exp(g_r)}
```

```math
E_{\mathrm{mix}} = \sum_{m \in \mathcal M} w_m E_m
```

In addition to direct fusion, `MultiModalFusion` uses a six-head sequence-to-structure cross-attention module, taking sequence embeddings as **queries**, graph-derived structure embeddings as **keys and values**, and combining the attention result with the direct weighted sum through a residual addition and layer normalization. This is how structural information can modify the initial fused representation through learned dependencies rather than only through the input scalar gates. The source implementation runs this cross-attention in FP32 on CUDA to avoid the recorded BF16 nonfinite-output problem, then casts the result back to the main network dtype.

### S19.3 Six stacked Viroformer blocks and 12-head biased self-attention

The encoder applies six `ViroformerBlock` instances in sequence, updating nucleotide states at each depth. Their input/output hidden width is 768. Each block contains 12 self-attention heads of size 64, derived from (768/12). A learned linear projection generates query, key and value tensors for each head. At query nucleotide (i), key nucleotide (j) and head (h), the content score is:

```math
s_{ij}^{(h)} = \frac{(q_i^{(h)})^T k_j^{(h)}}{\sqrt{64}}
```

Two independent learned, head-specific terms are added to that content score. The first depends on the conservation scores at **both** participating nucleotides and is capped below at zero before applying `log1p`. The second scales the UFold-derived pair-contact entry:

```math
b_{ij,\mathrm{cons}}^{(h)} = \ln\!\left(1 + \max\!\bigl[\alpha_h(c_i+c_j),0\bigr]\right)
```

```math
b_{ij,\mathrm{struct}}^{(h)} = \beta_h P_{ij}
```

The attention logit is `s + b_cons + b_struct`. Invalid key positions are masked before softmax; each head's normalized attention weights mix its value vectors, and the 12 head outputs are concatenated into width 768 and passed through a learned output projection. The trainable coefficients `alpha[h]` and `beta[h]` permit different heads to use conservation and structural evidence differently. These biases are **soft inductive preferences**, not hard constraints enforcing biological pairing or regulatory function. The implementation can use PyTorch scaled-dot-product attention for efficiency; explicitly requested attention weights use the diagnostic calculation path, so interpretability exports should document the sampled checkpoint and evaluation example.

### S19.4 Host cross-attention, normalization, regularization and feed-forward computation

After the self-attention residual update, each encoder block applies host cross-attention using the nucleotide representation as query and the broadcast 768-dimensional host representation as keys and values. The model then adds that output and applies layer normalization. Since all host positions share the same projected global vector, a near-uniform host attention distribution is architecturally expected; the heatmap alone does not indicate sequence-localized host binding.

A feed-forward subnetwork expands the per-position representation from 768 to 3,072 dimensions and projects it back to 768 with GELU, dropout and a residual-plus-normalization pathway. Source settings include dropout 0.10 and stochastic depth 0.10. In the implementation, stochastic depth can bypass an entire block's transform during training, while the encoder remains fully connected at evaluation. The host cross-attention operation is evaluated in FP32 under CUDA mixed-precision execution; the resulting update returns to the surrounding working dtype.

### S19.5 Refined pair state, gate and prediction outputs

The initial pair state is the UFold-derived contact matrix. The six-block encoder produces updated nucleotide states, from which the pair decoder evaluates pairs of nucleotide representations using a two-layer learned network. The resulting pair values pass through a sigmoid, are symmetrized by averaging with their transpose, and are thresholded. This yields a computationally refined, symmetric pair state rather than experimentally measured RNA contacts.

One refinement pass was configured for the production architecture. A structure gate computes a sigmoid of a learned function of the initial fused representation and a row-wise summary of the refined pair state. The final nucleotide state is `gate * initial_fused + (1 - gate) * encoder_updated`. This lets the network mix structural updates with original fused features. It does not establish that the refined structure is biologically correct.

The `MotifHeads` module maps the final states to motif-presence logits and probabilities, 87-category motif-type logits and probabilities, boundary/segmentation logits, conservation logits, predicted accessibility, host-related logits, and optional masked-language-model and motif-projection outputs. The motif-type head uses the available annotation classes; positions with no supported class use ignore index −1 during supervised loss construction. The downstream interval decoder operates on nucleotide prediction probabilities and applies frozen validation-defined thresholds before the motif catalog is assembled.

### S19.6 Learning objectives, optional branches, and what the experiment did not establish

The primary segmentation loss combines masked binary cross-entropy with sigmoid Dice loss, so learning is driven by both nucleotide-level classification and interval overlap. The implementation supports six groups of objectives under learned homoscedastic uncertainty weighting: conservation, motif segmentation, accessibility, host prediction, contrastive representation learning, and an auxiliary group. The 87-class motif-type cross-entropy and separate motif-presence BCE are added to the grouped objective. Label smoothing for motif type is 0.10. Accessibility and host losses are active only when valid corresponding target masks are present. Contrastive InfoNCE requires explicitly provided positive pairs, masked-language-model cross-entropy requires its input mask, denoising requires clean pair targets, and the KL invariance term requires a mutated latent. The source methodology reports that optional contrastive and KL-invariance terms were inactive in the standard training procedure unless the required inputs were supplied.

The presence of a branch in `virad.py` does not prove that it contributed to a particular recorded run. Confirm activation from training configuration or logs before claiming an ablation or a quantified improvement from that branch. Similarly, the staged 2+5+5+43 epoch curriculum is not a matched set of independently trained modality ablations, and attention entropy does not constitute causal evidence of RNA motif function.

### S19.7 Suggested internal architecture illustration

## S20. Complete Viroformer computational methodology


### S20.1 Why the model uses four biological inputs

Sequence alone can identify local nucleotide combinations but does not explicitly encode whether a site is evolutionarily constrained, whether distant nucleotides are predicted to contact one another, or which host-associated protein context was supplied. VIRAD therefore aligns **four sources to the same viral window**. The sequence tensor contains base identity; the phyloP vector represents evolutionary conservation scored from training-only family models; UFold supplies a binary predicted RNA contact matrix; ESM-2 supplies a global 1,280-dimensional protein-context vector. The source paper characterizes these as complementary representations for RNA-region discovery, **not** four independent experimental measurements.

The original tensor contract is:
| Quantity | Model-ready size | Interpretation and safeguard |
|---|---|---|
| Nucleotide window | `[B,5,512]` | One-hot A, U, G, C, N; accession and coordinates recorded |
| Conservation | `[B,512]` | Mapped phyloP conservation from frozen training-only model |
| RNA contacts | `[B,512,512]` | UFold output with finite, symmetric, binary semantics |
| Host-associated context | `[B,1280]` | Global ESM-2 embedding, not observed nucleotide contact |
| Position mask | `[B,512]` | Removes invalid positions from relevant computations |
| Motif presence | `[B,512]` | RNA-relevant label at each nucleotide |
| Motif type | `[B,512]` | One of 87 valid classes or ignore index −1 |
| Fused/encoder states | `[B,512,768]` | Shared learned representation |

**Algorithm 4** requires consistent sequence and conservation slices, the correct contact-cache provenance, finite numeric tensors, and aligned label coordinates. A mismatch rejects the example rather than silently assigning a conservation score or predicted contact to the wrong base. The 512-nt window with 256-nt stride is the *training* context; the later 1,024-nt genomic inference context is a separate processing stage.

### S20.2 Modality embedding before transformer attention

The **nucleotide projection** converts the five one-hot channels to the shared dimension 768; sinusoidal positional encoding supplies information about the relative location of a base inside its local window. The **conservation projection** converts each scalar score through a learned `1 -> 128 -> 768` network. The **structure projection** uses contact-aware message passing so that a nucleotide can aggregate features from other nucleotides connected to it by the binary UFold map, even if they are far apart in primary-sequence coordinates. The **host projection** transforms a single 1,280-component protein vector into 768 components and repeats this across all 512 nucleotide positions.

Each representation has compatible size `[B,512,768]` by the time learned fusion occurs. These projections **do not make the four inputs interchangeable**: sequence supplies base identity and order, conservation supplies evolutionary scores, the graph branch supplies predicted contact neighborhood, and the host branch supplies shared contextual conditioning.

### S20.3 Learned fusion and separate sequence-to-structure cross-attention

The paper specifies four learned modality gate logits, denoted `g_m`, and a normalized weight `w_m` for modality `m` in the set `M = {sequence, conservation, structure, host}`. Source **Equation 12** is:

```math
w_m = \frac{\exp(g_m)}{\sum_{r\in\mathcal M}\exp(g_r)}
```

Source **Equation 13** combines all four projected states:

```math
E_{\mathrm{mix}} = \sum_{m\in\mathcal M} w_m E_m
```

The softmax ensures the modality weights are nonnegative and sum to one. A larger weight indicates a stronger learned contribution in the trained fusion operation, **not a direct estimate of biological causal importance**. Direct fusion is supplemented by **six-head sequence-to-structure cross-attention**, where the sequence-derived states are queries and the structure-derived states supply keys and values. The resulting structural-context output is combined with the direct fusion through residual addition and normalization. This is important because the sum alone would mix already projected embeddings without explicitly allowing a sequence position to retrieve information from other positions in its structural context.

### S20.4 The exact content, conservation, and structural terms inside each Viroformer block

The model takes the fused `[B,512,768]` tensor through **six sequential transformer blocks**. With **12 self-attention heads** and total width 768, each head operates on width `d_h = 768/12 = 64`. The model learns query, key and value projections for each position and each head.

At nucleotide positions `i` and `j`, source **Equation 14** gives the ordinary scaled content-attention score:

```math
s_{ij}^{(h)} = \frac{\bigl(q_i^{(h)}\bigr)^T k_j^{(h)}}{\sqrt{64}}
```

This term is a learned similarity between position `i`'s query and position `j`'s key. It allows nucleotide positions far apart on the linear genome to influence one another. The originality of the VIRAD attention specification is that the raw content score is **not the only term**: learned biological biases enter before softmax.

For source **Equation 15**, the conservation-related bias uses the conservation values `c_i` and `c_j` of both nucleotides and a head-specific learned coefficient (alpha_h):

```math
b_{ij,\mathrm{cons}}^{(h)}
= \ln\!\left(1 + \max\!\left[\alpha_h(c_i+c_j),0\right]\right)
```

For source **Equation 16**, a head-specific coefficient `beta_h` scales the binary UFold contact `P_ij`:

```math
b_{ij,\mathrm{struct}}^{(h)} = \beta_h P_{ij}
```

The combined source **Equation 17** is:

```math
a_{ij}^{(h)}=s_{ij}^{(h)}
+b_{ij,\mathrm{cons}}^{(h)}
+b_{ij,\mathrm{struct}}^{(h)}
```

The model then masks invalid keys, normalizes the logits with softmax over eligible key positions, and uses the resulting weights to aggregate value states. A high structural bias means a predicted contact *may* be emphasized, not that the algorithm proves the contact occurs in infected cells. Each of the 12 heads has its own learned (alpha_h) and `beta_h`, so heads can differ in how much sequence similarity, evolutionary conservation, and predicted pairing contribute. Crucially, the conservation term depends on **both positions**, whereas the structure term describes the **particular pair** of positions.

### S20.5 Full order of computation within one of the six blocks

Source **Algorithm 5** describes each Viroformer block in the following order. First, the nucleotide states are normalized. Second, the block produces the 12 attention heads' query, key, and value tensors and adds the conservation and predicted-contact biases to their similarity scores. Third, the masked attention output returns a new `[B,512,768]` context representation and is combined with the block input through a residual update. Fourth, host cross-attention conditions the nucleotide states on the projected global host vector. Fifth, a positionwise feed-forward network expands width **768 to 3,072**, applies **GELU**, and projects back to **768**, with residual and normalization pathways.

Dropout was **0.10**, and stochastic depth was **0.10** in the reported production architecture. The same sequence of operations is repeated six times, allowing progressively transformed nucleotide states to depend on different positions and biological signals. Six blocks are **six successive transformations**, not six independent models or six separate viral families. The 12 heads within each block are learned parallel attention components, not 12 distinct biological tests.

### S20.6 Iterative pair-state refinement and structure gate

The first pair state comes from UFold's binary contact matrix. After the transformer has updated nucleotide representations, the dense pair decoder produces pairwise predictions from contextual states. The predicted scores are symmetrized so that the contact from `i` to `j` agrees with the contact from `j` to `i`. The study documentation states that the final model used **one structural-refinement pass**, rather than an unspecified number of iterative cycles.

A learned sigmoid structure gate controls the degree to which the encoder's structurally updated representation replaces the original fused representation. Let (E_{mathrm{mix}}) be the original fused state, `Z` the encoder-updated state, and (gamma) the gate. The implemented combination is:

```math
Z_{\mathrm{final}} = \gamma\odot E_{\mathrm{mix}}+(1-\gamma)\odot Z
```

The structure gate protects against **assuming predicted UFold contacts must dominate every position**. It does not make the refined map experimentally verified. The output projection heads subsequently produce per-nucleotide motif probability and boundary/segmentation values, available motif-type categories, and auxiliary targets used in training.

### S20.7 Multitask prediction targets and conditional learning losses

The framework's primary prediction problem is nucleotide-level motif localization. The source methodology distinguishes **motif-presence supervision** from a more specific **87-category motif-type task**. Invalid or unannotated type positions use ignore index −1, preventing missing labels from being interpreted as a real motif class.

The original methodological framework combines a positionwise **masked binary cross-entropy** objective with **Dice loss** for motif segmentation. BCE penalizes inaccurate nucleotide classification, whereas Dice loss rewards overlap of predicted and known motif regions, which is particularly relevant for sparse regions. The corresponding source equations are:

```math
\mathcal L_{\mathrm{Dice}}
=1-\frac{2\sum_i p_i y_i}
{\max\!\left(\sum_i p_i+\sum_i y_i,10^{-8}\right)}
```

```math
\mathcal L_{\mathrm{motif}}=\mathcal L_{\mathrm{Dice}}+\mathcal L_{\mathrm{BCE}}
```

The conservation auxiliary target is based on the threshold (c_i>0.7), with masked BCE for relevant positions. Motif-type categorical cross-entropy uses **label smoothing 0.10**, applied only where a type label exists. When corresponding labels are actually supplied, the framework can also optimize accessibility error, host-label prediction, masked-language-model reconstruction, denoising, and contrastive representations. However, optional tasks do not become active simply because a head exists in the model source code.

The source methodological specification uses a **six-group learned uncertainty weighting** for available tasks rather than arbitrarily adding all auxiliary losses at fixed weights:

```math
\widetilde{\mathcal L}_k=
\tfrac12\exp(-2s_k)\mathcal L_k+s_k
```

```math
\mathcal L_{\mathrm{group}}=
\sum_{k\in\mathcal A}\widetilde{\mathcal L}_k
```

Here, (mathcal{A}) is the set of active, supervised task groups and `s_k` is the learned task log-uncertainty parameter. Motif-type and motif-presence losses are added to that grouped objective; the optional KL-invariance term requires an explicitly supplied mutated-latent representation. The source methodology states that positive-pair contrastive and KL invariance terms were not active in the ordinary run when their required data were absent. **This distinction matters for reproducibility:** a described potential loss is not automatically proof that it influenced the reported fitted network.

### S20.8 Exact staged training, numerical stability, and model selection

The study's four-stage training specification defines a four-stage curriculum that progressively introduces biological information.

| Stage | Upper limit | Information available | Starting learning rate |
|---|---:|---|---:|
| 1 | 2 epochs | Nucleotide sequence | 1 × 10⁻⁴ |
| 2 | 5 epochs | Sequence, conservation, structure | 1 × 10⁻⁴ |
| 3 | 5 epochs | Above plus host context | 5 × 10⁻⁵ |
| 4 | 43 epochs | Full multimodal and multitask optimization | 1 × 10⁻⁵ |

The maximum configured training duration was 55 epochs, while the completed reported run reached **global epoch 51** and **95,064 optimizer steps**. Its AdamW settings were `beta_1 = 0.9`, `beta_2 = 0.999`, weight decay **0.01**, gradient clipping at norm **1.0**, linear warmup, and cosine learning-rate decay toward **10⁻⁶**. A micro-batch of **4** with **8 accumulated micro-batches** formed effective batch size **32**. The exponential-moving-average teacher used decay **0.999**, meaning its parameters were updated using 99.9% previous teacher and 0.1% current student values.

Training ran on the reported NVIDIA L40 46-GB configuration with BF16 mixed precision and TF32/SDPA support. The paper documents numerically unstable outputs in the sequence-to-structure and host cross-attention during stage four; those operations were therefore computed locally in **FP32** and cast back to the surrounding working type. This was a numerical-stability change, not a change to the six-block architecture.

**Checkpoint selection and early stopping used different criteria.** Strict improvement in validation loss selected the saved best checkpoint. Stage-four early stopping monitored validation F1 with patience **8** and minimum global epoch **20**. Final-state files were saved each epoch, while periodic checkpoints were kept every five global epochs, with recent checkpoints retained. Training examples supplied gradients, validation data controlled selection, and untouched test data were not used for optimization or hyperparameter choice.

### S20.9 Architectural interpretation and empirical boundary

The 12-head attention mechanism and six encoder blocks are computational architectural components. Learned attention weights describe model-internal weighting of nucleotide positions and do not experimentally establish RNA base pairing, biological function, or host-protein contacts.

A genuine claim that conservation, structure, or host context *causes* an improvement in prediction requires matched retraining comparisons and held-out scores, not merely stage-wise addition of inputs during one curriculum. Likewise, successful sequence-level Cas13d guide matching on unseen accessions does not demonstrate guide cleavage, antiviral treatment efficacy, or generalization to wholly unrepresented viral families. These restrictions are retained in this supplementary archive.
