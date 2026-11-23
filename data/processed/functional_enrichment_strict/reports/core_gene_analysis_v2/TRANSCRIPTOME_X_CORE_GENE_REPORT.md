# TRANSCRIPTOME-X

# Core gene analysis - Ensembl-based integration

## 1. Methodology

All significant GO and Reactome enrichment terms were used.
Genes were extracted from the readable geneID field.
Enrichment symbols were mapped to Ensembl stable identifiers.
Only unambiguous one-to-one symbol-Ensembl-Entrez mappings were retained.
Differential expression results were integrated using Ensembl stable IDs.
Enrichment and DE gene symbols were preserved in separate columns.

## 2. Annotation versions

AnnotationDbi: 1.74.0
org.Hs.eg.db: 3.23.1

## 3. Summary

       comparison direction core_genes genes_with_ensembl genes_matched_to_DE
 siNRP1_3_vs_siNT        UP          0                  0                   0
 siNRP1_3_vs_siNT      DOWN        744                703                 703
 siNRP1_5_vs_siNT        UP        861                809                 809
 siNRP1_5_vs_siNT      DOWN        156                143                 143
 significant_DE_genes direction_consistent_genes
                    0                          0
                  703                        703
                  809                        809
                  143                        143

## 4. Mapping status

    mapping_status symbol_count
  AMBIGUOUS_SYMBOL          100
 STRICT_ONE_TO_ONE         1557

## 5. DE integration status

              de_mapping_status gene_count
             MATCHED_BY_ENSEMBL       1655
 NO_UNAMBIGUOUS_ENSEMBL_MAPPING        106

## 6. Top 20 core genes

### siNRP1_3_vs_siNT - UP

No core genes identified.

### siNRP1_3_vs_siNT - DOWN

 enrichment_symbol de_gene_name      ensembl_id unique_terms           methods      logFC
              BMP4         BMP4 ENSG00000125378           62 BP;CC;MF;REACTOME -1.0747844
               NOG          NOG ENSG00000183691           62             BP;MF -0.6833364
            CTNNB1       CTNNB1 ENSG00000168036           58             BP;CC -1.6638434
              DAG1         DAG1 ENSG00000173402           46 BP;CC;MF;REACTOME -0.9668756
             TGFB2        TGFB2 ENSG00000092969           46    BP;CC;REACTOME -0.9879095
            TGFBR3       TGFBR3 ENSG00000069702           42          BP;CC;MF -1.9007557
              NRP1         NRP1 ENSG00000099250           39          BP;CC;MF -0.6928770
              FBN1         FBN1 ENSG00000166147           37 BP;CC;MF;REACTOME -0.9198911
              CDH5         CDH5 ENSG00000179776           33          BP;CC;MF -1.7181682
              PTEN         <NA>            <NA>           33             BP;CC         NA
             THBS1        THBS1 ENSG00000137801           32 BP;CC;MF;REACTOME -1.9130429
             SMAD6        SMAD6 ENSG00000137834           32                BP -1.3493933
              HES1         HES1 ENSG00000114315           32                BP -1.3044943
             PPARG        PPARG ENSG00000132170           31             BP;CC -0.6263427
               SKI          SKI ENSG00000157933           30                BP -1.0588117
              NRP2         NRP2 ENSG00000118257           30          BP;CC;MF -1.6534550
             ITGA2        ITGA2 ENSG00000164171           30 BP;CC;MF;REACTOME -0.7073159
               WT1          WT1 ENSG00000184937           29                BP -1.3645540
              CCN3         CCN3 ENSG00000136999           28          BP;CC;MF -0.9850759
              TBX1         TBX1 ENSG00000184058           28                BP -0.8442807
          FDR de_direction direction_consistent              de_mapping_status
 7.130926e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 3.024785e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 6.610548e-05         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 3.376873e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 4.455604e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 2.319161e-04         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 2.478353e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 6.007932e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 4.882182e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
 9.850147e-05         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 4.105804e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 1.228872e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 1.309479e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 7.247148e-04         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 1.103260e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 3.001419e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 4.632826e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 1.203308e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 1.784529e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL

### siNRP1_5_vs_siNT - UP

 enrichment_symbol de_gene_name      ensembl_id unique_terms           methods     logFC
              EDN1         EDN1 ENSG00000078401          110             BP;MF 0.8823842
              BMP4         BMP4 ENSG00000125378           85    BP;MF;REACTOME 0.7200604
              CAV1         CAV1 ENSG00000105974           84 BP;CC;MF;REACTOME 0.6660968
             TGFB2        TGFB2 ENSG00000092969           83 BP;CC;MF;REACTOME 1.0660136
             PDGFB        PDGFB ENSG00000100311           73 BP;CC;MF;REACTOME 2.8069988
              CAV3         CAV3 ENSG00000182533           73 BP;CC;MF;REACTOME 1.9727500
               KDR          KDR ENSG00000128052           70 BP;CC;MF;REACTOME 2.0710847
               JUP          JUP ENSG00000173801           70    BP;CC;REACTOME 0.8105077
              TPM1         TPM1 ENSG00000140416           67 BP;CC;MF;REACTOME 0.5932634
               ENG          ENG ENSG00000106991           66          BP;CC;MF 1.2101445
             PPARG        PPARG ENSG00000132170           65             BP;MF 1.2663597
            PRKACA         <NA>            <NA>           65 BP;CC;MF;REACTOME        NA
             UBA52        UBA52 ENSG00000221983           62 BP;CC;MF;REACTOME 0.5470508
             MEF2C        MEF2C ENSG00000081189           62       BP;REACTOME 1.0937635
            SEMA5A       SEMA5A ENSG00000112902           60    BP;MF;REACTOME 2.0435373
             RPL39        RPL39 ENSG00000198918           59 BP;CC;MF;REACTOME 0.6928513
               MDK          MDK ENSG00000110492           59 BP;CC;MF;REACTOME 0.9133485
            SLC9A1       SLC9A1 ENSG00000090020           58 BP;CC;MF;REACTOME 0.6443364
             RPS19        RPS19 ENSG00000105372           58 BP;CC;MF;REACTOME 0.5298886
              NOS3         NOS3 ENSG00000164867           56 BP;CC;MF;REACTOME 0.6598733
          FDR de_direction direction_consistent              de_mapping_status
 4.954609e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 3.560301e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 4.085667e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 1.559315e-03           UP                 TRUE             MATCHED_BY_ENSEMBL
 7.880129e-05           UP                 TRUE             MATCHED_BY_ENSEMBL
 4.027009e-03           UP                 TRUE             MATCHED_BY_ENSEMBL
 8.565578e-05           UP                 TRUE             MATCHED_BY_ENSEMBL
 2.197782e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 3.560301e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 9.615606e-04           UP                 TRUE             MATCHED_BY_ENSEMBL
 9.115681e-05           UP                 TRUE             MATCHED_BY_ENSEMBL
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
 2.105644e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 3.553099e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 1.834058e-03           UP                 TRUE             MATCHED_BY_ENSEMBL
 3.967448e-03           UP                 TRUE             MATCHED_BY_ENSEMBL
 1.423886e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 7.483334e-03           UP                 TRUE             MATCHED_BY_ENSEMBL
 3.999571e-02           UP                 TRUE             MATCHED_BY_ENSEMBL
 2.702151e-02           UP                 TRUE             MATCHED_BY_ENSEMBL

### siNRP1_5_vs_siNT - DOWN

 enrichment_symbol de_gene_name      ensembl_id unique_terms        methods      logFC
          HLA-DRB5     HLA-DRB5 ENSG00000198502           13 BP;CC;REACTOME -2.2009611
          HLA-DPA1         <NA>            <NA>           13 BP;CC;REACTOME         NA
          HLA-DPB1         <NA>            <NA>           13 BP;CC;REACTOME         NA
           HLA-DRA         <NA>            <NA>           13 BP;CC;REACTOME         NA
          HLA-DRB1         <NA>            <NA>           13 BP;CC;REACTOME         NA
           HLA-DMA         <NA>            <NA>           12          BP;CC         NA
           HLA-DMB         <NA>            <NA>           12          BP;CC         NA
             HLA-F         <NA>            <NA>            6 BP;CC;REACTOME         NA
               HFE          HFE ENSG00000010704            4             BP -3.2131742
              CTSS         CTSS ENSG00000163131            4             BP -1.9023239
             MFSD6        MFSD6 ENSG00000151690            4             BP -1.4759925
             LNPEP        LNPEP ENSG00000113441            4             BP -1.3979293
            TAPBPL       TAPBPL ENSG00000139192            4             BP -1.2063449
              CTSL         CTSL ENSG00000135047            4             BP -0.4841283
              CD44         CD44 ENSG00000026508            3    BP;REACTOME -0.9732231
              CD74         CD74 ENSG00000019582            3          BP;CC -1.0921939
             ULBP2        ULBP2 ENSG00000131015            2             BP -2.6373770
            ACTL6A       ACTL6A ENSG00000136518            2             BP -2.2567003
            SLAIN1       SLAIN1 ENSG00000139737            2             CC -3.4104541
              IL18         IL18 ENSG00000150782            2             BP -2.0574271
          FDR de_direction direction_consistent              de_mapping_status
 2.010717e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
           NA         <NA>                FALSE NO_UNAMBIGUOUS_ENSEMBL_MAPPING
 8.213386e-06         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 2.727514e-05         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 2.468518e-04         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 1.518189e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 2.872470e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 4.294615e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 3.551760e-03         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 1.141629e-02         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 7.681506e-06         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 9.318494e-06         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 2.388343e-05         DOWN                 TRUE             MATCHED_BY_ENSEMBL
 3.514684e-05         DOWN                 TRUE             MATCHED_BY_ENSEMBL

## 7. Important limitations

- Gene recurrence is not an independent measure of biological importance.
- GO terms are hierarchical and may be redundant.
- Reactome pathways may share many genes.
- Ambiguous symbol mappings are not assigned arbitrarily.
- Genes without an unambiguous mapping remain in recurrence tables but are not integrated with DE.
- The same dataset is used for differential expression and enrichment.
- Recurrence does not establish causality or pathway activation.
- Expression direction is checked against the original DE results.

## 8. Output files

- gene_pathway_membership_ensembl.csv
- gene_symbol_mapping_audit.csv
- gene_symbol_mapping_strict.csv
- core_gene_frequency_by_method.csv
- core_gene_frequency_overall.csv
- core_genes_with_differential_expression.csv
- top_20_core_genes.csv
- mapping_status_summary.csv
- DE_integration_status_summary.csv
- core_gene_summary.csv
