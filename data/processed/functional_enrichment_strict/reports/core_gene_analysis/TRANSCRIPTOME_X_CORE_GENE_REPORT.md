# TRANSCRIPTOME-X

# Core gene analysis

## 1. Methodology

All significant GO and Reactome enrichment terms were used.
Genes were extracted from the readable geneID field.
Gene recurrence was calculated using unique enrichment term IDs.
Core genes were integrated with differential expression results.

## 2. Summary

       comparison direction core_genes significant_DE_genes direction_consistent_genes
 siNRP1_3_vs_siNT        UP          0                    0                          0
 siNRP1_3_vs_siNT      DOWN        744                  742                        742
 siNRP1_5_vs_siNT        UP        861                  860                        860
 siNRP1_5_vs_siNT      DOWN        156                  156                        156

## 3. Top 20 core genes

### siNRP1_3_vs_siNT - UP

No core genes identified.

### siNRP1_3_vs_siNT - DOWN

 gene_name unique_terms           methods      logFC          FDR de_direction
      BMP4           62 BP;CC;MF;REACTOME -1.0747844 7.130926e-03         DOWN
       NOG           62             BP;MF -0.6833364 3.024785e-02         DOWN
    CTNNB1           58             BP;CC -1.6638434 6.610548e-05         DOWN
      DAG1           46 BP;CC;MF;REACTOME -0.9668756 3.376873e-03         DOWN
     TGFB2           46    BP;CC;REACTOME -0.9879095 4.455604e-03         DOWN
    TGFBR3           42          BP;CC;MF -1.9007557 2.319161e-04         DOWN
      NRP1           39          BP;CC;MF -0.6928770 2.478353e-02         DOWN
      FBN1           37 BP;CC;MF;REACTOME -0.9198911 6.007932e-03         DOWN
      PTEN           33             BP;CC -0.5784064 3.350612e-02         DOWN
      CDH5           33          BP;CC;MF -1.7181682 4.882182e-02         DOWN
     THBS1           32 BP;CC;MF;REACTOME -1.9130429 9.850147e-05         DOWN
     SMAD6           32                BP -1.3493933 4.105804e-03         DOWN
      HES1           32                BP -1.3044943 1.228872e-02         DOWN
     PPARG           31             BP;CC -0.6263427 1.309479e-02         DOWN
       SKI           30                BP -1.0588117 7.247148e-04         DOWN
      NRP2           30          BP;CC;MF -1.6534550 1.103260e-03         DOWN
     ITGA2           30 BP;CC;MF;REACTOME -0.7073159 3.001419e-02         DOWN
       WT1           29                BP -1.3645540 4.632826e-03         DOWN
      CCN3           28          BP;CC;MF -0.9850759 1.203308e-02         DOWN
      TBX1           28                BP -0.8442807 1.784529e-02         DOWN
 direction_consistent
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE

### siNRP1_5_vs_siNT - UP

 gene_name unique_terms           methods     logFC          FDR de_direction
      EDN1          110             BP;MF 0.8823842 4.954609e-02           UP
      BMP4           85    BP;MF;REACTOME 0.7200604 3.560301e-02           UP
      CAV1           84 BP;CC;MF;REACTOME 0.6660968 4.085667e-02           UP
     TGFB2           83 BP;CC;MF;REACTOME 1.0660136 1.559315e-03           UP
     PDGFB           73 BP;CC;MF;REACTOME 2.8069988 7.880129e-05           UP
      CAV3           73 BP;CC;MF;REACTOME 1.9727500 4.027009e-03           UP
       KDR           70 BP;CC;MF;REACTOME 2.0710847 8.565578e-05           UP
       JUP           70    BP;CC;REACTOME 0.8105077 2.197782e-02           UP
      TPM1           67 BP;CC;MF;REACTOME 0.5932634 3.560301e-02           UP
       ENG           66          BP;CC;MF 1.2101445 9.615606e-04           UP
     PPARG           65             BP;MF 1.2663597 9.115681e-05           UP
    PRKACA           65 BP;CC;MF;REACTOME 0.5090313 2.105644e-02           UP
     UBA52           62 BP;CC;MF;REACTOME 0.5470508 2.105644e-02           UP
     MEF2C           62       BP;REACTOME 1.0937635 3.553099e-02           UP
    SEMA5A           60    BP;MF;REACTOME 2.0435373 1.834058e-03           UP
     RPL39           59 BP;CC;MF;REACTOME 0.6928513 3.967448e-03           UP
       MDK           59 BP;CC;MF;REACTOME 0.9133485 1.423886e-02           UP
    SLC9A1           58 BP;CC;MF;REACTOME 0.6443364 7.483334e-03           UP
     RPS19           58 BP;CC;MF;REACTOME 0.5298886 3.999571e-02           UP
      NOS3           56 BP;CC;MF;REACTOME 0.6598733 2.702151e-02           UP
 direction_consistent
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE

### siNRP1_5_vs_siNT - DOWN

 gene_name unique_terms        methods      logFC          FDR de_direction
   HLA-DRA           13 BP;CC;REACTOME -2.2127136 1.037523e-03         DOWN
  HLA-DPB1           13 BP;CC;REACTOME -2.0608376 2.308310e-03         DOWN
  HLA-DPA1           13 BP;CC;REACTOME -2.0408979 3.943277e-03         DOWN
  HLA-DRB1           13 BP;CC;REACTOME -1.9135809 1.132135e-02         DOWN
  HLA-DRB5           13 BP;CC;REACTOME -2.2009611 2.010717e-02         DOWN
   HLA-DMA           12          BP;CC -1.8010297 9.573789e-04         DOWN
   HLA-DMB           12          BP;CC -2.1113836 1.112846e-03         DOWN
     HLA-F            6 BP;CC;REACTOME -1.6961280 1.598168e-03         DOWN
       HFE            4             BP -3.2131742 8.213386e-06         DOWN
      CTSS            4             BP -1.9023239 2.727514e-05         DOWN
     MFSD6            4             BP -1.4759925 2.468518e-04         DOWN
     LNPEP            4             BP -1.3979293 1.518189e-03         DOWN
    TAPBPL            4             BP -1.2063449 2.872470e-03         DOWN
      CTSL            4             BP -0.4841283 4.294615e-02         DOWN
      CD44            3    BP;REACTOME -0.9732231 3.551760e-03         DOWN
      CD74            3          BP;CC -1.0921939 1.141629e-02         DOWN
     ULBP2            2             BP -2.6373770 7.681506e-06         DOWN
    ACTL6A            2             BP -2.2567003 9.318494e-06         DOWN
    SLAIN1            2             CC -3.4104541 2.388343e-05         DOWN
      IL18            2             BP -2.0574271 3.514684e-05         DOWN
 direction_consistent
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE
                 TRUE

## 4. Important limitations

- Gene recurrence is not an independent measure of biological importance.
- GO terms are often hierarchical and redundant.
- Reactome pathways may share many genes.
- A gene occurring in many related terms may reflect annotation structure.
- Readable gene symbols can be ambiguous; ambiguous symbol-to-Ensembl mappings are reported separately.
- The same dataset is used for differential expression and enrichment.
- Recurrence does not establish causality or pathway activation.
- Expression direction is checked against the original differential expression results.

## 5. Output files

- gene_pathway_membership.csv
- core_gene_frequency_by_method.csv
- core_gene_frequency_overall.csv
- core_genes_with_differential_expression.csv
- top_20_core_genes.csv
- ambiguous_gene_symbols.csv
- core_gene_summary.csv
