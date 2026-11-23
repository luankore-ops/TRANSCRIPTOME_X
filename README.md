# TRANSCRIPTOME-X
## Integrated Transcriptomic, Functional and Network-Based Analysis of NRP1 Knockdown in Claudin-Low Breast Cancer Cells

**Project type:** Bioinformatics · RNA-seq · Differential Gene Expression · Functional Enrichment · Systems Biology · Network Analysis  
**Dataset:** GEO GSE266566  
**Biological model:** MDA-MB-231 claudin-low breast cancer cells  
**Primary research focus:** Transcriptional and functional changes associated with Neuropilin-1 (NRP1) knockdown  
**Analysis environment:** R 4.6.2 · Bioconductor 3.23  
**Project status:** Completed

---

# 1. Project Overview

TRANSCRIPTOME-X is an integrated bioinformatics project designed to investigate transcriptional changes associated with Neuropilin-1 (NRP1) knockdown in MDA-MB-231 claudin-low breast cancer cells.

The project combines differential gene expression analysis, normalization sensitivity assessment, functional enrichment, pathway-level characterization, gene-module analysis, biological network reconstruction, candidate gene prioritization and integrated biological interpretation.

Rather than focusing exclusively on individual differentially expressed genes, the analytical framework progressively integrates multiple layers of molecular evidence to identify transcriptional patterns and candidate genes supported by complementary analytical approaches.

Two independent NRP1-targeting knockdown conditions are investigated and compared against a non-targeting control. Particular attention is given to the consistency of transcriptional responses across the two knockdown conditions, the functional context of the identified genes and their relationships within biological networks.

The final objective is to develop a reproducible and interpretable analytical workflow capable of connecting differential expression results with higher-level functional and network-based evidence.

# 2. Scientific Background

Neuropilin-1 (NRP1) is a multifunctional transmembrane co-receptor involved in several biological processes, including growth factor signalling, cellular communication, adhesion and tissue organization.

Its biological relevance extends to cancer, where alterations in NRP1-associated signalling may be related to changes in cellular behaviour and interactions with the surrounding microenvironment.

Claudin-low breast cancer represents a molecular subtype characterized by distinctive transcriptional features, including alterations in cell adhesion, epithelial differentiation and extracellular matrix-associated processes.

The MDA-MB-231 cell line is frequently used as an experimental model for investigating aggressive breast cancer phenotypes and associated molecular mechanisms.

Within this context, transcriptomic analysis following NRP1 knockdown provides an opportunity to characterize molecular expression patterns associated with reduced NRP1 expression and to investigate the functional processes represented by the resulting transcriptional changes.

The present study adopts a comparative and integrative approach, examining two NRP1-targeting conditions independently before combining their results into a broader analytical framework.

# 3. Research Objectives

The main objectives of TRANSCRIPTOME-X are:

1. Characterize transcriptional changes associated with NRP1 knockdown through differential gene expression analysis.
2. Evaluate the robustness of differential expression results across alternative normalization strategies.
3. Identify enriched biological processes, molecular functions, cellular components and signalling pathways.
4. Compare functional signatures between the two knockdown conditions.
5. Investigate the relationship between differentially expressed genes and biologically relevant gene modules.
6. Reconstruct integrated gene-pathway-module networks to characterize molecular connectivity.
7. Identify candidate genes supported by multiple complementary evidence layers.
8. Develop an integrated biological interpretation of the observed transcriptional patterns.
9. Generate publication-oriented visualizations to communicate the principal findings.

# 4. Dataset and Experimental Design

## 4.1 Data source

The analysis uses the publicly available RNA-seq dataset:

**GEO accession:** GSE266566

**Publication:**

Lynam LR et al. (2026). *Bulk RNA sequencing dataset of Claudin-low breast cancer cell lines with Neuropilin-1 knockdown.*

Scientific Data, 13, 20.

DOI: https://doi.org/10.1038/s41597-025-06332-7

The dataset consists of bulk RNA-sequencing measurements obtained from MDA-MB-231 claudin-low breast cancer cells under different experimental conditions.

## 4.2 Experimental groups

The dataset includes nine biological samples distributed across three experimental groups.

| Experimental condition | GEO samples | Biological replicates |
|---|---|---:|
| Non-targeting control (siNT) | GSM8250651–GSM8250653 | 3 |
| NRP1 knockdown (siNRP1.3) | GSM8250645–GSM8250647 | 3 |
| NRP1 knockdown (siNRP1.5) | GSM8250648–GSM8250650 | 3 |
| **Total** | **9 samples** | **9** |

The two knockdown conditions represent different NRP1-targeting reagents.

The following comparisons were considered:

- siNRP1.3 versus siNT.
- siNRP1.5 versus siNT.

These comparisons were analysed separately to preserve the distinction between the two experimental conditions.

## 4.3 Expression data

The original expression data consist of RSEM expected counts.

These values are fractional rather than necessarily integer-valued.

Because effective gene lengths could not be reconstructed from the available input data, the analysis retained the supplied expected counts and used library-size normalization approaches implemented in edgeR.

The original expression values were not artificially rounded to integers, and no unavailable gene-length information was inferred or fabricated.

This methodological choice and its implications are considered when interpreting the results.

# 5. Computational Environment

The project was developed using R and Bioconductor.

| Component | Version |
|---|---|
| R | 4.6.2 |
| Bioconductor | 3.23 |
| edgeR | 4.10.7 |
| clusterProfiler | 4.20.3 |
| ReactomePA | 1.56.0 |
| AnnotationDbi | 1.74.0 |
| org.Hs.eg.db | 3.23.1 |

The analytical workflow uses R-based statistical and bioinformatics tools for expression analysis, gene annotation, functional enrichment and biological interpretation.

# 6. Analytical Workflow

TRANSCRIPTOME-X follows a multistage analytical strategy.

The workflow progressively transforms gene-level expression measurements into integrated biological evidence.

## 6.1 Differential Gene Expression Analysis

Differential expression analysis was performed independently for the two NRP1 knockdown conditions.

The primary comparisons were:

- siNRP1.3 versus siNT.
- siNRP1.5 versus siNT.

The analysis used edgeR with TMM normalization.

Differential expression results were evaluated using log2 fold changes and false discovery rate-adjusted p-values.

Genes with FDR below 0.05 were considered statistically significant for the principal differential expression analyses.

### Main results

| Comparison | Significant genes |
|---|---:|
| siNRP1.3 vs siNT | 2,690 |
| siNRP1.5 vs siNT | 3,420 |

A comparison of the two differential expression signatures identified:

- 737 genes significant in both conditions.
- 533 genes with concordant expression directions.
- 204 genes with discordant expression directions.

These findings indicate that the two knockdown conditions share a subset of transcriptional responses while also exhibiting condition-specific and directionally different expression patterns.

## 6.2 Normalization Sensitivity Analysis

To evaluate the dependence of differential expression results on normalization choices, alternative normalization approaches were investigated.

The strategies included:

- TMM.
- TMMwsp.
- RLE.
- Upper-quartile normalization.

The resulting differential expression estimates were compared using log-fold-change correlations and overlap-based measures, including the Jaccard index.

The analysis showed high agreement between normalization strategies, supporting the relative robustness of the observed transcriptional patterns to the normalization method.

This assessment does not eliminate limitations associated with the original expression data, but provides an additional measure of analytical consistency.

## 6.3 Functional Enrichment Analysis

Functional enrichment analyses were performed to investigate the biological processes and pathways represented by the differential expression signatures.

The analysis included functional annotation and pathway-oriented interpretation using available gene annotation resources and Bioconductor packages.

The principal biological themes included:

- Extracellular matrix organization and remodelling.
- Glycosylation-associated processes.
- Cell adhesion and cell-substrate interactions.
- Growth factor-associated signalling.
- Translation-related processes.
- Cytoskeletal and contractile-associated components.

The enrichment results were interpreted as evidence of functional representation within the analysed gene signatures.

Enrichment does not, by itself, demonstrate pathway activation or inhibition.

## 6.4 Integrated Pathway Analysis

Pathway-level results from the two knockdown conditions were compared to characterize shared and condition-specific functional signatures.

The pathway analysis identified 282 distinct pathways.

Their distribution was:

| Pathway category | Number |
|---|---:|
| Shared | 12 |
| Specific to siNRP1.3 | 39 |
| Specific to siNRP1.5 | 231 |
| **Total** | **282** |

The pathway prioritization framework further classified pathways according to descriptive evidence-support categories:

| Evidence category | Number |
|---|---:|
| High | 44 |
| Medium | 198 |
| Exploratory | 40 |

The prioritization score integrated multiple descriptive components, including gene support, significant gene fraction, directional concordance and module support.

These categories represent relative analytical evidence and should not be interpreted as independent statistical significance levels.

## 6.5 Gene-Module Analysis

A gene-module analytical layer was introduced to investigate the relationships between candidate genes and biologically relevant functional modules.

Four principal modules were considered:

- ECM_REMODELING.
- TGFB_BMP_SIGNALING.
- CELL_ADHESION_CYTOSKELETON.
- GROWTH_FACTOR_RECEPTOR_SIGNALING.

The integrated analysis identified 155 unique genes across 159 gene-module assignments.

Among these genes:

- 40 were significant in both knockdown comparisons.
- 27 showed concordant expression directions.
- 13 showed discordant expression directions.

Several genes were associated with multiple modules, suggesting that they may represent connections between different functional contexts.

Examples include CCN2, CCN3 and THBS1.

## 6.6 Integrated Gene-Pathway-Module Network

A multilayer biological network was reconstructed by integrating:

- Gene nodes.
- Pathway nodes.
- Functional module nodes.
- Gene-pathway relationships.
- Gene-module relationships.
- Pathway-module relationships.

The integrated network contained:

| Network component | Number |
|---|---:|
| Nodes | 2,199 |
| Edges | 15,556 |
| Connected components | 5 |
| Louvain communities | 16 |
| Articulation points | 72 |

Network topology was investigated using connectivity and centrality-related measures, including degree and betweenness centrality.

The analysis highlighted genes with potentially relevant positions within the integrated network.

Among the genes with high betweenness centrality were:

- TGFB2.
- THBS1.
- BMP4.
- CAV1.
- CDH1.
- ANXA1.
- LYN.
- NRP1.
- CTNNB1.
- IL6R.

These measurements describe network topology and do not independently establish biological causality or functional importance.

## 6.7 Network-Based Candidate Prioritization

A network-based prioritization strategy was used to rank candidate genes according to their relative importance within the reconstructed biological framework.

A total of 1,913 candidate genes were evaluated.

The resulting prioritization categories were:

| Category | Number of genes |
|---|---:|
| High priority | 75 |
| Medium priority | 886 |
| Exploratory | 952 |
| **Total** | **1,913** |

The ranking was subsequently integrated with differential expression, pathway and gene-module evidence.

## 6.8 Integrated Candidate Gene Validation

The integrated validation stage combined four complementary evidence layers:

1. Network-based prioritization.
2. Differential expression in both knockdown comparisons.
3. Pathway-level support.
4. Gene-module support.

A total of 1,913 candidate genes were evaluated.

The final integrated classification was:

| Evidence classification | Number | Percentage |
|---|---:|---:|
| Robust cross-layer | 23 | 1.20% |
| Multi-layer supported | 56 | 2.93% |
| Transcriptomically supported | 1,810 | 94.62% |
| Network priority only | 0 | 0% |
| Exploratory | 24 | 1.25% |
| **Total** | **1,913** | **100%** |

The final shortlist contained 79 candidate genes.

The 23 robust cross-layer candidates were characterized by statistical significance in both knockdown comparisons, concordant expression direction and additional pathway and gene-module support.

Importantly, the term *robust cross-layer* refers to convergence across the analytical evidence layers available in this dataset. It does not indicate independent experimental validation.

# 7. Principal Results

## 7.1 Final candidate prioritization

The integrated prioritization analysis identified the following leading candidates:

| Rank | Gene | Integrated score |
|---:|---|---:|
| 1 | THBS1 | 97.41 |
| 2 | CCN3 | 97.08 |
| 3 | NRP1 | 94.63 |
| 4 | TGFBR3 | 94.51 |
| 5 | VEGFC | 94.04 |
| 6 | PIK3CA | 94.02 |
| 7 | FGF1 | 93.76 |
| 8 | RELN | 93.39 |
| 9 | MET | 92.57 |
| 10 | PAK2 | 92.50 |

These genes represent candidates supported by the integrated analytical framework and provide a focused set for further biological investigation.

## 7.2 Concordant transcriptional patterns

Among the 79 shortlisted candidates, the expression patterns were distributed as follows:

| Expression pattern | Number | Percentage |
|---|---:|---:|
| Downregulated in both conditions | 15 | 18.99% |
| Upregulated in both conditions | 8 | 10.13% |
| Opposite directions | 13 | 16.46% |
| Significant only in siNRP1.3 | 16 | 20.25% |
| Significant only in siNRP1.5 | 27 | 34.18% |
| **Total** | **79** | **100%** |

All 23 robust cross-layer candidates showed concordant expression directions:

- 15 were downregulated in both knockdown conditions.
- 8 were upregulated in both knockdown conditions.

This convergence provides a particularly relevant subset for subsequent investigation.

## 7.3 Principal robust candidates

The 23 robust cross-layer candidates were:

**Downregulated in both conditions**

THBS1, CCN3, NRP1, TGFBR3, PIK3CA, COL4A2, FBLN1, COL5A1, PLOD2, LUM, MATN2, HTRA1, FRAS1, TIMP4.

**Upregulated in both conditions**

VEGFC, FGF1, RELN, MET, CCN2, MMP1, ADAMTS6, TAGLN.

These candidates are associated with several biological contexts, particularly extracellular matrix organization, cell signalling, cellular adhesion and structural organization.

Their expression patterns should be interpreted within the limitations of the experimental model and the available transcriptomic data.

## 7.4 Biological interpretation

The integrated results suggest that NRP1 knockdown is associated with transcriptional changes involving multiple interconnected biological processes.

Three broad functional themes emerge:

**1. Extracellular matrix organization and remodelling**

Several prioritized genes, including THBS1, CCN3, CCN2, COL4A2, COL5A1, FBLN1 and PLOD2, are associated with extracellular matrix-related processes.

Their representation across differential expression, functional annotation and network analyses supports the relevance of this biological context within the dataset.

**2. Growth factor and signalling-associated processes**

Genes such as NRP1, TGFBR3, VEGFC, FGF1 and MET are represented among the prioritized candidates.

Their presence suggests that growth factor-associated signalling constitutes another relevant functional dimension of the observed transcriptional response.

**3. Cell adhesion and structural organization**

The functional enrichment and network analyses also identify relationships involving cellular adhesion, cytoskeletal organization and cell-substrate interactions.

Together, these findings suggest that the transcriptional responses associated with NRP1 knockdown involve interconnected extracellular matrix, signalling and cellular organization-related processes.

These observations are hypothesis-generating and require further experimental investigation.

# 8. Final Visualizations

The final visualization stage generated seven publication-oriented figures.

Each figure was exported in both PNG and PDF formats, with PNG files produced at 300 dpi.

| Figure | Description |
|---|---|
| Figure 01 | Ranking of the 23 robust cross-layer candidate genes |
| Figure 02 | Expression heatmap of the 79 shortlisted candidates |
| Figure 03 | Integrated evidence matrix |
| Figure 04 | Robust gene-module network |
| Figure 05 | Distribution of integrated candidate classifications |
| Figure 06 | Comparison of log2 fold changes between knockdown conditions |
| Figure 07 | Top 15 integrated candidate genes |

## Output directories

**PNG figures**

`results/final_results_visualization/png/`

**PDF figures**

`results/final_results_visualization/pdf/`

**Audit and validation outputs**

`results/final_results_visualization/audit/`

The figures provide complementary representations of candidate prioritization, expression patterns, evidence integration and network relationships.

# 9. Reproducibility and Output Organization

The project is organized into sequential analytical stages, with intermediate and final results stored in dedicated directories.

The principal output directories include:

```text
TRANSCRIPTOME_X/
│
├── results/
│   ├── normalization_sensitivity/
│   ├── functional_enrichment_interpretation/
│   ├── pathway_concordance_specificity/
│   ├── gene_pathway_module_network/
│   ├── network_candidate_prioritization/
│   ├── integrated_candidate_gene_validation/
│   ├── final_biological_interpretation/
│   └── final_results_visualization/
│       ├── png/
│       ├── pdf/
│       └── audit/
│
└── scripts/
```

The directory structure above summarizes the principal analytical outputs and is intended as a guide; the actual project directory remains the authoritative reference for script and file organization.

Intermediate tables, reports, prioritization results and audit outputs are retained to facilitate inspection of the analytical decisions and the resulting evidence.

# 10. Methodological Limitations

Several limitations must be considered when interpreting the results.

## 10.1 Limited sample size

Each experimental group contains three biological replicates.

Although this design allows differential expression analysis, the limited number of samples reduces statistical power and may affect the stability of estimates, particularly for genes with modest expression changes.

## 10.2 Single cellular model

The analysis is based on the MDA-MB-231 cell line.

Consequently, the findings cannot automatically be generalized to other breast cancer cell lines, tumour subtypes or patient populations.

## 10.3 Interpretation of RSEM expected counts

The available expression measurements consist of fractional RSEM expected counts.

Effective gene lengths could not be reconstructed from the available inputs. The analysis therefore retained the supplied values and applied library-size normalization.

This is an important methodological constraint when interpreting the differential expression estimates.

## 10.4 Knockdown-condition comparability

The two NRP1-targeting conditions represent different reagents evaluated against the same control group.

Although concordant expression patterns provide useful evidence of consistency, the two comparisons should not be considered independent biological replications.

Differences between the knockdown conditions may reflect differences in knockdown efficiency, reagent-specific effects or other experimental factors.

## 10.5 Functional enrichment limitations

Enrichment analyses identify overrepresented biological annotations within a gene set.

They do not directly measure pathway activity, establish pathway activation or inhibition, or demonstrate a causal relationship between NRP1 and the identified biological processes.

## 10.6 Network interpretation

The reconstructed network integrates relationships derived from the available analytical evidence and functional annotations.

Network centrality and connectivity measures describe the topology of this representation.

They should not be interpreted as direct experimental measurements of molecular interactions or proof that a particular gene controls a biological process.

## 10.7 Candidate prioritization scores

The integrated scores are descriptive ranking measures.

They summarize the convergence of different evidence layers and are not equivalent to statistical significance, effect size, clinical relevance or experimentally demonstrated biological importance.

The final robust candidate classification represents internal analytical consistency rather than independent experimental validation.

# 11. Conclusions

TRANSCRIPTOME-X provides an integrated computational framework for investigating transcriptomic changes associated with NRP1 knockdown in MDA-MB-231 claudin-low breast cancer cells.

The analysis identified substantial differential expression signatures in both knockdown conditions and demonstrated that a subset of genes exhibits concordant transcriptional responses.

The integration of differential expression, normalization sensitivity, functional enrichment, pathway characterization, gene-module analysis and network topology enabled the identification of 79 shortlisted candidates, including 23 robust cross-layer candidates.

The principal biological themes emerging from the analysis involve extracellular matrix organization and remodelling, growth factor-associated signalling, cell adhesion and cellular structural organization.

Among the prioritized candidates, THBS1, CCN3, NRP1 and TGFBR3 showed particularly strong integrated evidence.

Overall, the project demonstrates the value of combining gene-level statistical analysis with functional and network-based approaches to obtain a more comprehensive interpretation of transcriptomic datasets.

The results provide a structured basis for generating biological hypotheses and selecting candidate genes for future experimental validation.

However, the findings remain computational and exploratory. Independent experimental studies, additional cellular models and complementary molecular assays would be required to establish the functional significance and potential mechanistic relationships of the identified candidates.

# 12. Future Perspectives

Potential future developments include:

- Experimental validation of selected candidate genes using qPCR or other molecular assays.
- Validation of protein-level changes through Western blotting or immunofluorescence.
- Investigation of extracellular matrix-related phenotypes.
- Functional characterization of selected growth factor-associated signalling pathways.
- Replication of the analysis in additional breast cancer cell lines.
- Integration with independent transcriptomic datasets.
- Comparison with publicly available clinical and tumour expression datasets.
- Investigation of whether the prioritized signatures are associated with clinically relevant breast cancer phenotypes.

These extensions would help determine the reproducibility, biological significance and potential translational relevance of the computational findings.

# 13. References

1. Lynam LR et al. (2026). *Bulk RNA sequencing dataset of Claudin-low breast cancer cell lines with Neuropilin-1 knockdown.* Scientific Data, 13, 20. https://doi.org/10.1038/s41597-025-06332-7

2. Robinson MD, McCarthy DJ, Smyth GK. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

3. McCarthy DJ, Chen Y, Smyth GK. (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Research*, 40(10), 4288–4297. https://doi.org/10.1093/nar/gks042

4. Wu T et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141

5. Ritchie ME et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47. https://doi.org/10.1093/nar/gkv007

---

**TRANSCRIPTOME-X — Integrated Bioinformatics and Systems Biology Analysis**

*An end-to-end computational investigation of NRP1-associated transcriptional responses in claudin-low breast cancer cells.*