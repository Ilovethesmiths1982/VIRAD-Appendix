# VIRAD Supplementary Figure Legends

This document provides supplementary figure descriptions and notes on the available model-derived visualization artifacts.

## Supplementary Figure S4. Attention heatmaps from the six Viroformer blocks

**[INSERT SUPPLEMENTARY FIGURE S4 IMAGE HERE , six aligned self-attention panels (a–f) plus cross-attention panel (g). Prefer an uncluttered 2 × 3 arrangement for (a–f) and a separate full-width (g), rather than seven tiny panels.]**

**Supplementary Figure S4. Self-attention and cross-attention patterns in VIRAD.** **(a–f)** Heatmaps of model-derived self-attention for Viroformer Blocks 1, 2, 3, 4, 5, and 6, respectively, across nucleotide positions in a 512-nt window. **(g)** Cross-attention patterns derived from the multimodal architecture. The panels visualize how attention weights are distributed within and between model representations and complement the quantitative attention-entropy and effective-position summaries shown in main Figure 3. These attention weights are **computational model diagnostics**, not direct experimental evidence of RNA pairing, motif function, or host–virus interaction. **Source-paper aggregation:** Original research Figure 7 describes mean self-attention heatmaps averaged across the 12 attention heads in each of the six Viroformer blocks. Original Figure 8 shows fusion and host cross-attention. **The local PNG exports have not been opened or verified**, so confirm that they correspond to these source-paper figures before asserting that the files themselves use the same aggregation. Source accession, genomic coordinates, checkpoint, evaluation partition, and color normalization should be checked in the export metadata.

### Exact local source files

| Panel | Filename reported in local PowerShell search |
|---|---|
| **S4a** | `Figure_10_Block_1_self_attention_heatmap.png` |
| **S4b** | `Figure_10_Block_2_self_attention_heatmap.png` |
| **S4c** | `Figure_10_Block_3_self_attention_heatmap.png` |
| **S4d** | `Figure_10_Block_4_self_attention_heatmap.png` |
| **S4e** | `Figure_10_Block_5_self_attention_heatmap.png` |
| **S4f** | `Figure_10_Block_6_self_attention_heatmap.png` |
| **S4g** | `Figure_11_cross_attention_patterns.png` |

**Local folder reported:** `C:\Users\Zhyrex\Downloads\New folder\`. These **local images have not been uploaded to this repository**. Do not replace them with illustrative or invented heatmaps. One candidate array file is `attention_matrices.npz` in the user's VIRAD figure directories; use it only after checking its tensor shapes, heads, coordinate provenance and checkpoint.

### Supplementary Results paragraph

The six Viroformer blocks exhibited distinct self-attention distributions (Supplementary Figure S4a–f), indicating differences in how successive encoder layers weighted nucleotide positions within the input window. Cross-attention patterns (Supplementary Figure S4g) provide an additional view of interactions between the model's learned representations. The quantitative analysis in the main article showed relatively diffuse self-attention across much of the 512-nt window and comparatively concentrated fusion attention; host cross-attention was near-uniform under the recorded settings. These findings describe the model's internal behavior and are not sufficient, on their own, to infer causal biological mechanisms.

## Interpretation and availability

**Main Figure 3** summarizes normalized attention entropy and effective attended positions across the six self-attention blocks, fusion attention and host cross-attention. **Supplementary Figure S4(a–g)** shows the detailed heatmaps. Using one grouped S4 caption prevents each block from counting as a separate main-paper display item.

The final supplementary PDF or Word attachment should place the actual graphic **immediately above its S4 caption**. Confirm the event's rules for supplementary figures and the readability of labels when printed.
