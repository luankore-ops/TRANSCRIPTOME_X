# TRANSCRIPTOME-X

**RNA-seq differential expression and functional enrichment analysis of NRP1 knockdown in MDA-MB-231 breast cancer cells.**

## Overview

TRANSCRIPTOME-X is a bioinformatics analysis of public RNA-seq dataset GSE266566. It investigates transcriptional changes associated with two independent NRP1 knockdown conditions (siNRP1.3 and siNRP1.5) compared with non-targeting control (siNT) in MDA-MB-231 cells.

The analysis includes differential expression, cross-knockdown agreement, functional enrichment, and gene-level annotation and expression-profile quality checks.

## Dataset and experimental design

- **Source:** [NCBI GEO – GSE266566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266566)
- **Cell line:** MDA-MB-231
- **Conditions:** siNT, siNRP1.3, siNRP1.5
- **Biological replicates:** 3 per condition (9 samples total)
- **Input:** RSEM expected counts

The input matrix contains 78,894 gene rows and 9 samples. Lowly expressed genes were filtered with edgeR::filterByExpr; 17,531 genes were retained. Fractional expected counts were retained without rounding.

## Methods

- TMM normalization
- Differential expression with edgeR quasi-likelihood methodology and robust dispersion estimation
- Separate contrasts: siNRP1_3_vs_siNT and siNRP1_5_vs_siNT
- Significance threshold: FDR < 0.05
- Functional enrichment using Gene Ontology (BP, CC, MF) and Reactome
- Cross-knockdown concordance analysis
- Gene identity/annotation and duplicate expression-profile audits

## Main results

### Differential expression

| Contrast | DE genes (FDR < 0.05) |
|---|---:|
| siNRP1.3 vs siNT | 2,690 |
| siNRP1.5 vs siNT | 3,420 |

NRP1 (ENSG00000099250.18) showed negative differential expression in both contrasts:

| Contrast | logFC | FDR |
|---|---:|---:|
| siNRP1.3 vs siNT | -0.692877 | 0.02478353 |
| siNRP1.5 vs siNT | -4.876953 | 3.023637e-07 |

These are model-based RNA-seq estimates, not direct measurements of protein-level knockdown efficiency.

### Cross-knockdown concordance

Among 17,531 genes tested in both contrasts, 737 were significant in both; 533 showed the same direction of change and 204 showed opposite directions. The concordant set is a prioritization result, not independent proof of reproducibility.

### Functional enrichment

The strict enrichment analysis identified different patterns across contrasts and directions. Notable results included extracellular-matrix organization and glycosylation among genes downregulated in siNRP1.3, and translation-related, focal-adhesion/cell-substrate-junction, and contractile/cytoskeletal annotations among genes upregulated in siNRP1.5.

Enrichment represents statistical over-representation of annotations. It does not by itself demonstrate pathway activation, inhibition, or causal mechanisms. Related GO and Reactome terms may overlap.

## Quality checks and interpretation

The project includes audits of gene identifiers, annotation mappings, cross-contrast direction consistency, and identical expression profiles. Identical profiles or ambiguous annotations are reported rather than automatically collapsed or removed.

Results are specific to this dataset and experimental model. Follow-up validation is required before drawing mechanistic or translational conclusions.

## Reproducibility

The analysis was performed in R 4.6.2 with Bioconductor 3.23. Key package versions include edgeR 4.10.7, clusterProfiler 4.20.3, ReactomePA 1.56.0, AnnotationDbi 1.74.0, and org.Hs.eg.db 3.23.1.

Raw data are not redistributed; retrieve the original dataset from GEO using the accession above. Analysis scripts, selected result tables, figures, and additional documentation will be added from the local project directory.

## Author

**Luan Kore**