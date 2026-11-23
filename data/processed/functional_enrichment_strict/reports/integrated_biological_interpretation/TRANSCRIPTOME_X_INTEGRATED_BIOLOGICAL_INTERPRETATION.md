# TRANSCRIPTOME-X

# Integrated biological interpretation

## 1. Study overview

Dataset: GSE266566.

Cell model: MDA-MB-231.

Experimental conditions:
- siNT: non-targeting control
- siNRP1_3: NRP1 knockdown condition 3
- siNRP1_5: NRP1 knockdown condition 5

Comparisons:
- siNRP1_3 vs siNT
- siNRP1_5 vs siNT

Differential expression was assessed using edgeR quasi-likelihood testing.
Functional enrichment uses strict one-to-one Ensembl-Entrez mapping.

## 2. Differential expression summary

       comparison tested_genes significant_genes up_regulated down_regulated
 siNRP1_3_vs_siNT        17531              2690         1165           1525
 siNRP1_5_vs_siNT        17531              3420         1647           1773

## 3. NRP1 expression

       comparison            gene_id gene_name     logFC          FDR fold_change
 siNRP1_3_vs_siNT ENSG00000099250.18      NRP1 -0.692877 2.478353e-02  0.61861899
 siNRP1_5_vs_siNT ENSG00000099250.18      NRP1 -4.876953 3.023637e-07  0.03403227

## 4. Cross-comparison enrichment overlap

 direction   method terms_siNRP1_3 terms_siNRP1_5 shared_terms only_siNRP1_3 only_siNRP1_5
        UP       BP              0            269            0             0           269
        UP       MF              0             20            0             0            20
        UP       CC              0             45            0             0            45
        UP REACTOME              0             60            0             0            60
      DOWN       BP            110             10            0           110            10
      DOWN       MF             14              0            0            14             0
      DOWN       CC             29              4            0            29             4
      DOWN REACTOME              4              1            0             4             1
 jaccard
       0
       0
       0
       0
       0
       0
       0
       0

## 5. Most significant enriched terms

### siNRP1_3_vs_siNT - UP - BP

No significantly enriched terms were identified.

### siNRP1_3_vs_siNT - UP - MF

No significantly enriched terms were identified.

### siNRP1_3_vs_siNT - UP - CC

No significantly enriched terms were identified.

### siNRP1_3_vs_siNT - UP - REACTOME

No significantly enriched terms were identified.

### siNRP1_3_vs_siNT - DOWN - BP

         ID                                          Description     p.adjust Count
 GO:0009100                       glycoprotein metabolic process 2.486362e-05    62
 GO:0009101                    glycoprotein biosynthetic process 3.007279e-05    54
 GO:0006493                       protein O-linked glycosylation 7.820776e-05    33
 GO:0030198                    extracellular matrix organization 4.999465e-04    41
 GO:0043062                 extracellular structure organization 4.999465e-04    41
 GO:0045229        external encapsulating structure organization 4.999465e-04    41
 GO:0003013                           circulatory system process 1.152108e-03    62
 GO:0001503                                         ossification 1.311977e-03    49
 GO:0003205                          cardiac chamber development 1.311977e-03    22
 GO:0008285 negative regulation of cell population proliferation 1.311977e-03    72

### siNRP1_3_vs_siNT - DOWN - MF

         ID                          Description     p.adjust Count
 GO:0098631      cell adhesion mediator activity 0.0002761077    23
 GO:0019838                growth factor binding 0.0013369827    27
 GO:0019955                     cytokine binding 0.0025962265    23
 GO:0008194     UDP-glycosyltransferase activity 0.0025962265    26
 GO:0016757         glycosyltransferase activity 0.0040949168    40
 GO:0016758          hexosyltransferase activity 0.0084696635    31
 GO:0005509                  calcium ion binding 0.0097321486    67
 GO:0098632 cell-cell adhesion mediator activity 0.0122760260    14
 GO:0005518                     collagen binding 0.0122760260    15
 GO:0008201                      heparin binding 0.0122760260    22

### siNRP1_3_vs_siNT - DOWN - CC

         ID                      Description     p.adjust Count
 GO:0045178               basal part of cell 0.0002697310    47
 GO:0031012             extracellular matrix 0.0006049368    56
 GO:0030312 external encapsulating structure 0.0006049368    56
 GO:0009925            basal plasma membrane 0.0007961416    42
 GO:0005788      endoplasmic reticulum lumen 0.0013470180    47
 GO:0098984         neuron to neuron synapse 0.0013754075    53
 GO:0016323      basolateral plasma membrane 0.0037054388    36
 GO:0005796                      Golgi lumen 0.0052401681    16
 GO:0005911               cell-cell junction 0.0052401681    65
 GO:0070160                   tight junction 0.0052447824    25

### siNRP1_3_vs_siNT - DOWN - REACTOME

            ID                        Description     p.adjust Count
 R-HSA-1474244  Extracellular matrix organization 0.0002823214    49
 R-HSA-5173105             O-linked glycosylation 0.0175544297    22
  R-HSA-216083 Integrin cell surface interactions 0.0175544297    18
 R-HSA-3000178                  ECM proteoglycans 0.0333675476    16

### siNRP1_5_vs_siNT - UP - BP

         ID                                  Description     p.adjust Count
 GO:0003012                        muscle system process 1.307482e-07    67
 GO:0006941                  striated muscle contraction 1.155658e-06    41
 GO:0061061                 muscle structure development 1.155658e-06    82
 GO:0006936                           muscle contraction 1.697765e-06    53
 GO:0050678  regulation of epithelial cell proliferation 1.697765e-06    48
 GO:0001936 regulation of endothelial cell proliferation 2.584786e-06    30
 GO:0001944                      vasculature development 6.959458e-06    86
 GO:0006935                                   chemotaxis 1.196869e-05    58
 GO:0042330                                        taxis 1.196869e-05    58
 GO:1901342        regulation of vasculature development 1.196869e-05    49

### siNRP1_5_vs_siNT - UP - MF

         ID                           Description     p.adjust Count
 GO:0003735    structural constituent of ribosome 1.923705e-07    44
 GO:0030546 signaling receptor activator activity 4.712273e-07    61
 GO:0030545 signaling receptor regulator activity 4.712273e-07    64
 GO:0048018              receptor ligand activity 6.439689e-07    59
 GO:0048306     calcium-dependent protein binding 5.937634e-06    20
 GO:0008083                growth factor activity 7.937851e-04    23
 GO:0005125                     cytokine activity 3.237434e-03    25
 GO:0019838                 growth factor binding 4.998811e-03    24
 GO:0003779                         actin binding 1.830097e-02    56
 GO:0004857             enzyme inhibitor activity 2.537703e-02    53

### siNRP1_5_vs_siNT - UP - CC

         ID                       Description     p.adjust Count
 GO:0030055           cell-substrate junction 1.817306e-10    86
 GO:0005925                    focal adhesion 5.401921e-10    83
 GO:0022626                cytosolic ribosome 1.715034e-08    36
 GO:0022627 cytosolic small ribosomal subunit 3.118602e-06    17
 GO:0044391                 ribosomal subunit 3.118602e-06    42
 GO:0043292          contractile muscle fiber 1.810426e-05    39
 GO:0030016                         myofibril 1.810426e-05    38
 GO:0030017                         sarcomere 2.478839e-05    35
 GO:0005911                cell-cell junction 7.503032e-05    70
 GO:0098803         respiratory chain complex 1.147368e-04    27

### siNRP1_5_vs_siNT - UP - REACTOME

            ID                                                   Description     p.adjust
  R-HSA-156842                             Eukaryotic Translation Elongation 4.068428e-07
 R-HSA-2408557                                      Selenocysteine synthesis 4.068428e-07
   R-HSA-72764                            Eukaryotic Translation Termination 4.068428e-07
  R-HSA-156902                                      Peptide chain elongation 4.068428e-07
  R-HSA-192823                                        Viral mRNA Translation 4.068428e-07
 R-HSA-9954714 PELO:HBS1L and ABCE1 dissociate a ribosome on a non-stop mRNA 4.068428e-07
   R-HSA-72689                      Formation of a pool of free 40S subunits 4.068428e-07
   R-HSA-72613                             Eukaryotic Translation Initiation 6.532737e-07
   R-HSA-72737                          Cap-dependent Translation Initiation 6.532737e-07
 R-HSA-9754678               SARS-CoV-2 modulates host translation machinery 1.358528e-06
 Count
    29
    29
    29
    28
    28
    28
    30
    33
    33
    20

### siNRP1_5_vs_siNT - DOWN - BP

         ID                                                      Description    p.adjust Count
 GO:0019882                              antigen processing and presentation 0.007653089    24
 GO:0048002           antigen processing and presentation of peptide antigen 0.027244047    18
 GO:0002399                            MHC class II protein complex assembly 0.028384014     7
 GO:0002503       peptide antigen assembly with MHC class II protein complex 0.028384014     7
 GO:0002478 antigen processing and presentation of exogenous peptide antigen 0.040322229    12
 GO:1903037                       regulation of leukocyte cell-cell adhesion 0.040322229    48
 GO:0019884         antigen processing and presentation of exogenous antigen 0.040322229    12
 GO:0007159                                     leukocyte cell-cell adhesion 0.040322229    52
 GO:0002396                                     MHC protein complex assembly 0.045560174     8
 GO:0002501                peptide antigen assembly with MHC protein complex 0.045560174     8

### siNRP1_5_vs_siNT - DOWN - MF

No significantly enriched terms were identified.

### siNRP1_5_vs_siNT - DOWN - CC

         ID                  Description    p.adjust Count
 GO:0042613 MHC class II protein complex 0.002387875     8
 GO:0005874                  microtubule 0.014045243    73
 GO:0042611          MHC protein complex 0.014045243     9
 GO:1990752              microtubule end 0.021623684    12

### siNRP1_5_vs_siNT - DOWN - REACTOME

           ID                Description  p.adjust Count
 R-HSA-877300 Interferon gamma signaling 0.0380719    22

## 6. Interpretation guidelines

- Enrichment indicates overrepresentation, not pathway activation or inhibition.
- Differential expression does not establish a causal mechanism.
- Results describe transcriptional associations in the studied cell model.
- GO terms may be redundant or biologically overlapping.
- The two knockdown conditions are compared separately against siNT.
- Their results do not constitute a direct statistical comparison between siNRP1_3 and siNRP1_5.
- The two comparisons use the same dataset and are not independent biological replications.
- Enrichment results depend on annotation coverage and the selected gene background.

## 7. Limitations

- RNA-seq data are RSEM expected counts, including fractional values.
- No RSEM inferential replicates were available to model read-assignment uncertainty.
- Strict one-to-one mapping excludes unmapped and ambiguous gene identifiers.
- Functional enrichment is based on genes passing the selected significance threshold.

## 8. Output files

- integrated_enrichment_results.csv
- top_enriched_terms.csv
- cross_comparison_enrichment_overlap.csv
- NRP1_expression_summary.csv
- differential_expression_summary.csv

This report is a computational summary. Biological conclusions require critical inspection of the enriched terms, their contributing genes, and the relevant literature.
