# TRANSCRIPTOME-X

**Integrated RNA-seq analysis and candidate prioritization of NRP1 knockdown in MDA-MB-231 claudin-low breast cancer cells.**

## Project overview

TRANSCRIPTOME-X is a computational reanalysis of the MDA-MB-231 subset of the public RNA-sequencing dataset [GEO GSE266566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266566). It investigates transcriptional changes associated with two independent NRP1-targeting siRNA conditions (siNRP1#3 and siNRP1#5) compared with a non-targeting control (siNT).

The project combines differential expression, cross-knockdown comparison, normalization sensitivity, functional enrichment, gene/network/module integration, candidate prioritization and biological hypothesis generation. The final computational workflow is documented through Script 51.

## Dataset and experimental design

- **Source:** [NCBI GEO – GSE266566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266566)
- **Cell line analyzed:** MDA-MB-231
- **Conditions:** siNT, siNRP1#3, siNRP1#5
- **Biological replicates:** 3 per condition (9 samples total)
- **Input:** RSEM gene-level expected counts
- **Input matrix:** 78,894 gene rows; 17,531 retained after expression filtering

The parent GEO series includes multiple claudin-low breast cancer cell lines and perturbation modalities. This repository's conclusions are restricted to the selected MDA-MB-231 siRNA subset.

## Methods

- Expression filtering with `edgeR::filterByExpr`
- TMM normalization for the primary differential-expression workflow
- Differential expression using edgeR quasi-likelihood methodology
- Separate contrasts: siNRP1#3 vs siNT and siNRP1#5 vs siNT
- Significance threshold: FDR < 0.05
- Normalization sensitivity comparison: TMM, TMMwsp, RLE and upper-quartile
- Functional enrichment using Gene Ontology and Reactome annotations
- Cross-knockdown concordance, annotation and expression-profile audits
- Integrated network, pathway and module evidence for candidate prioritization

## Main results

### Differential expression

| Contrast | DE genes (FDR < 0.05) |
|---|---:|
| siNRP1#3 vs siNT | 2,690 |
| siNRP1#5 vs siNT | 3,420 |

NRP1 showed negative differential expression in both contrasts. These are model-based RNA-seq estimates, not direct measurements of protein-level knockdown efficiency.

### Final integrated candidate assessment

| Metric | Result |
|---|---:|
| Genes evaluated in integrated candidate ranking | 1,913 |
| Final candidate shortlist | 79 |
| Candidate–pathway associations | 1,810 |
| Distinct pathway annotations | 247 |
| Candidate–module associations | 83 |
| Distinct modules | 4 |
| Research hypotheses generated | 21 |

Top-ranked candidates were THBS1, CCN3, NRP1, TGFBR3 and VEGFC. NRP1 ranked third, with 69 pathway associations and one module association (GROWTH_FACTOR_RECEPTOR_SIGNALING).

Pathway and module association counts reflect annotation coverage in the integrated evidence tables. They do not establish pathway activation, inhibition, direct molecular interaction or causality. Later evidence layers also reuse information generated in earlier analysis stages, so convergence should not be interpreted as fully independent validation.

## Final documentation

The final scientific dossier summarizes the study design, methods, quantitative results, candidate ranking, NRP1-focused interpretation, research hypotheses, limitations and proposed experimental validation:

- [Final scientific dossier (Markdown)](publication/reports/FINAL_SCIENTIFIC_DOSSIER.md)

The curated publication snapshot includes selected analysis scripts, result tables and figures under `publication/R`, `publication/tables` and `publication/figures`. It is not a byte-for-byte mirror of the local working directory.

## Reproducibility and data availability

The analysis was performed in R 4.6.2 with Bioconductor 3.23. Key package versions recorded in the project include edgeR 4.10.7, clusterProfiler 4.20.3, ReactomePA 1.56.0, AnnotationDbi 1.74.0 and org.Hs.eg.db 3.23.1.

Raw input data are not redistributed. Retrieve the original dataset from GEO using accession GSE266566. Local environments, caches, credentials, raw/external data and temporary files are excluded from this repository.

## Limitations

The analyzed subset has three biological replicates per condition and represents one cell model and one perturbation modality. Normalization robustness is not independent biological replication. Candidate ranking and pathway/module annotations support hypothesis generation, not mechanistic proof. Experimental validation is required before causal or translational conclusions can be drawn.

## Author

**Luan Kore**
