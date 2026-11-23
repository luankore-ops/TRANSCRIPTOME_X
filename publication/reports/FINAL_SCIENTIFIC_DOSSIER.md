# TRANSCRIPTOME-X — Final scientific dossier

## Project overview

TRANSCRIPTOME-X is a computational reanalysis of the MDA-MB-231 subset of GEO **GSE266566**, investigating transcriptional changes associated with two independent NRP1-targeting siRNA conditions (**siNRP1#3** and **siNRP1#5**) compared with a non-targeting control (**siNT**).

The analyzed subset contains three biological replicates per condition (nine samples total). The project combines RNA-seq quality and annotation checks, differential expression, cross-knockdown comparison, normalization sensitivity, functional annotation, network/module integration, candidate prioritization and biological hypothesis generation.

**Scope:** conclusions apply to this dataset and the analyzed MDA-MB-231 siRNA subset. They should not be generalized automatically to other cell lines, perturbation modalities or clinical settings.

## Dataset and primary analysis

- **Accession:** [GEO GSE266566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266566)
- **Organism / model:** *Homo sapiens* / MDA-MB-231
- **Groups:** siNT, siNRP1#3, siNRP1#5
- **Replicates:** 3 biological replicates per group
- **Input:** gene-level RSEM expected counts
- **Input features:** 78,894 gene rows
- **Genes retained after expression filtering:** 17,531
- **Differential expression:** edgeR quasi-likelihood workflow, with TMM as the primary normalization
- **Contrasts:** siNRP1#3 vs siNT; siNRP1#5 vs siNT
- **Significance threshold:** FDR < 0.05

### Differential-expression results

| Contrast | Significant genes (FDR < 0.05) |
|---|---:|
| siNRP1#3 vs siNT | 2,690 |
| siNRP1#5 vs siNT | 3,420 |

These are model-based RNA-seq estimates. They are not direct measurements of protein-level knockdown efficiency.

## Integrated analysis and candidate prioritization

The downstream workflow integrated differential-expression evidence with functional annotations, gene-network evidence and module assignments. The final candidate workflow evaluated **1,913 genes** and produced a **79-gene shortlist**.

The final evidence synthesis (Scripts 49–51) reported:

| Metric | Result |
|---|---:|
| Candidates assessed | 79 |
| Candidate–pathway associations | 1,810 |
| Distinct pathway annotations | 247 |
| Candidate–module associations | 83 |
| Distinct modules | 4 |
| Final candidate panel | 20 |
| Research hypotheses generated | 21 |
| Script 49 consistency checks | 6/6 passed |
| Script 50 consistency checks | 8/8 passed |
| Script 51 consistency checks | 7/7 passed |

### Top 20 prioritized candidates

| Rank | Gene | Pathway associations | Module associations |
|---:|---|---:|---:|
| 1 | THBS1 | 65 | 2 |
| 2 | CCN3 | 45 | 2 |
| 3 | NRP1 | 69 | 1 |
| 4 | TGFBR3 | 59 | 1 |
| 5 | VEGFC | 41 | 1 |
| 6 | PIK3CA | 39 | 1 |
| 7 | FGF1 | 39 | 1 |
| 8 | RELN | 30 | 1 |
| 9 | MET | 25 | 1 |
| 10 | PAK2 | 23 | 1 |
| 11 | COL4A2 | 19 | 1 |
| 12 | FBLN1 | 18 | 1 |
| 13 | CCN2 | 9 | 3 |
| 14 | COL5A1 | 16 | 1 |
| 15 | PLOD2 | 8 | 1 |
| 16 | LUM | 8 | 1 |
| 17 | MMP1 | 7 | 1 |
| 18 | MATN2 | 6 | 1 |
| 19 | ADAMTS6 | 6 | 1 |
| 20 | HTRA1 | 4 | 1 |

Pathway/module association counts describe annotation coverage in the assembled evidence tables. They are not measures of effect size, statistical significance, pathway activity or causal importance.

## NRP1-focused interpretation

NRP1 ranked **third** in the final candidate assessment, with **69 pathway associations** and one module association: **GROWTH_FACTOR_RECEPTOR_SIGNALING**.

This result makes growth-factor/receptor-associated signalling a reasonable area for follow-up. It does not demonstrate that this module is activated, inhibited or causally altered by NRP1 knockdown.

The leading candidates also include THBS1, CCN3, TGFBR3, VEGFC and FGF1, together with extracellular-matrix-associated genes such as COL4A2, COL5A1, FBLN1, PLOD2, LUM, MMP1, MATN2 and ADAMTS6. Collectively, these annotations motivate investigation of cell–matrix/extracellular-matrix biology and growth-factor/receptor-associated signalling. They do not establish a single unified mechanism.

## Normalization sensitivity

The project compared TMM, TMMwsp, RLE and upper-quartile normalization. The sensitivity analysis found substantial agreement between the resulting significant-gene sets and strong directional agreement among shared significant genes. This supports computational robustness to the tested normalization choices; it is not independent biological replication.

## Research hypotheses

1. **NRP1-associated transcriptional response:** the two NRP1-targeting siRNA contrasts are associated with broad transcriptional changes in the analyzed MDA-MB-231 model.
2. **Extracellular-matrix and cell–matrix response:** prioritized genes motivate follow-up of extracellular-matrix organization and cell–matrix-associated processes.
3. **Growth-factor and receptor-associated signalling:** NRP1, TGFBR3, VEGFC, FGF1, PIK3CA and MET, together with the NRP1-associated module annotation, motivate targeted investigation of receptor/growth-factor signalling.
4. **Shared versus perturbation-specific effects:** comparing the two siRNA contrasts can help distinguish shared signals from sequence-specific or perturbation-specific effects.

These are testable hypotheses generated from computational evidence, not experimentally established mechanisms.

## Limitations

- The analyzed subset contains only three biological replicates per condition.
- This project focuses on one cell model and one perturbation modality.
- Two siRNA sequences provide a useful comparison, but sequence-specific off-target effects remain possible.
- The main downstream analysis uses the public gene-level count data; it does not constitute a new independent sequencing experiment.
- Normalization sensitivity is not independent replication.
- Integrated evidence layers are not fully independent: later candidate assessments reuse evidence and indicators generated in earlier analysis stages.
- Pathway and module associations are annotation-based; they do not establish activation, inhibition, direct interaction or causality.
- The integrated candidate ranking is a prioritization framework, not a probability or formal validation score.
- Protein-structure visualizations are supplementary exploratory assets, not functional validation. The PyMOL rendering branch was frozen and is not presented as a validated result.

## Recommended experimental validation

1. Confirm NRP1 knockdown at RNA and protein level in the same cell model and experimental context.
2. Validate selected transcript changes by RT-qPCR using independent biological replicates.
3. Select a compact candidate panel using the underlying logFC, FDR, direction concordance, biological rationale and assay feasibility—not pathway-association counts alone.
4. Assess selected proteins with suitable assays (e.g. immunoblotting, flow cytometry or immunostaining).
5. Use an orthogonal NRP1 perturbation or rescue experiment to strengthen attribution to NRP1.
6. Perform phenotype-specific assays only after candidate-level validation (e.g. migration, invasion, adhesion or matrix-remodelling assays).

## Overall conclusion

TRANSCRIPTOME-X establishes an integrated computational workflow for prioritizing genes associated with NRP1-targeting siRNA conditions in the MDA-MB-231 subset of GSE266566. The analysis identifies substantial differential-expression signals in both contrasts and generates a 79-gene shortlist from a broader candidate space of 1,913 genes. THBS1, CCN3, NRP1 and TGFBR3 emerge among the leading candidates; NRP1 ranks third and is annotated to the GROWTH_FACTOR_RECEPTOR_SIGNALING module.

The strongest defensible conclusion is that the project produces an annotation-supported candidate panel and testable biological hypotheses. It does **not** demonstrate a causal NRP1 signalling mechanism, pathway activation or therapeutic effect. Independent experimental validation is required.

## Data and references

- Lynam L, Rockstroh A, Lehman M, Nelson C, Volpert M, Hollier B. *Bulk RNA sequencing dataset of Claudin-low breast cancer cell lines with Neuropilin-1 knockdown*. Scientific Data (2025). DOI: [10.1038/s41597-025-06332-7](https://doi.org/10.1038/s41597-025-06332-7).
- GEO accession: [GSE266566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266566).
- edgeR: Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*. 2010;26(1):139–140. DOI: [10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616).

---

*Prepared from the TRANSCRIPTOME-X project execution summaries and the curated public analysis outputs. Numerical values should be interpreted within the scope and limitations described above.*
