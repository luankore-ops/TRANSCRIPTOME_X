# TRANSCRIPTOME-X - Gene Identity and Annotation Audit

## Metodo

- Chiave primaria: Ensembl stable ID
- Rimozione suffisso versione Ensembl
- Annotazione tramite AnnotationDbi e org.Hs.eg.db
- Soglia DE: FDR < 0.05

## Risultati principali

- Geni concordanti analizzati: 533
- Geni che superano la validazione DE: 533
- Simboli corrispondenti all'annotazione: 487
- Simboli con mismatch: 9
- Geni senza annotazione: 36
- Simboli associati a piu Ensembl ID: 3

## Caso CRLF2

- NA | NA | NA
- ENSG00000205755 | CRLF2 | CRLF2
- ENSG00000292363 | CRLF2 | CRLF2
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA
- NA | NA | NA

## Controllo NRP1

- siNRP1_3_vs_siNT: logFC = -0.69288; FDR = 0.024784; direzione = DOWN
- siNRP1_5_vs_siNT: logFC = -4.877; FDR = 3.0236e-07; direzione = DOWN

## Interpretazione

Un simbolo genico associato a piu Ensembl ID non implica
necessariamente un errore: possono esistere identificativi
distinti che condividono un simbolo o rappresentano annotazioni
diverse. Gli identificativi non vengono collassati.

Le differenze tra gene_name e org.Hs.eg.db devono essere
interpretate considerando le versioni e le fonti di annotazione.

La corrispondenza con org.Hs.eg.db e un controllo di annotazione,
non una validazione sperimentale dell'identita biologica.

La validazione dei 533 geni utilizza gli stessi risultati DE
di R24 e non costituisce una replica sperimentale.
