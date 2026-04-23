# Gene Set Definitions

Standard gene sets used in Disease × AF module score analysis. All validated against Seurat RNA assay rownames from GSE166992.

## Inflammatory (31 genes)

Cytokines/Chemokines: IL6, IL1B, TNF, CXCL8, CCL2, CCL3, CCL4, CXCL1, CXCL2, CXCL3, CXCL10
NF-κB targets: NFKBIA, NFKB1, REL, RELA
IFN response: ISG15, ISG20, MX1, MX2, OAS1, IFIT1, IFIT2, IFIT3, IRF1, STAT1
Stress/acute phase: FOS, JUN, EGR1, ATF3, SOCS3, TNFAIP3

## Mitochondrial (15 genes)

13 MT-encoded proteins: MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, MT-CO1, MT-CO2, MT-CO3, MT-ATP6, MT-ATP8, MT-CYB
MT rRNA: MT-RNR1, MT-RNR2

## OXPHOS (33 genes)

Complex I: NDUFA1, NDUFA2, NDUFA4, NDUFB8, NDUFB10, NDUFS1, NDUFS2, NDUFS3, NDUFS7, NDUFS8
Complex II: SDHA, SDHB, SDHC, SDHD
Complex III: UQCRC1, UQCRC2, UQCRB, UQCRFS1
Complex IV: COX4I1, COX5A, COX5B, COX6A1, COX6B1, COX7A2, COX7C
Complex V: ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5MC1, ATP5MC2, ATP5MC3, ATP5MG

## Apoptosis (23 genes)

Pro-apoptotic: BAX, BAK1, BID, BIM*, BAD, PUMA*, NOXA1, CASP3, CASP7, CASP8, CASP9, CYCS, APAF1, DIABLO
Anti-apoptotic: BCL2, BCL2L1, MCL1
Death receptors: FAS, FASLG, TNFRSF10A, TNFRSF10B
Executioner: DFFA, DFFB

*BIM and PUMA are common aliases; may be missing from Seurat object (actual gene symbols: BCL2L11, BBC3)

## Type I IFN Response (29 genes)

ISGs: ISG15, ISG20, MX1, MX2, OAS1, OAS2, OAS3, IFIT1, IFIT2, IFIT3, IFIT5, IFITM1, IFITM2, IFITM3, IFI6, IFI27, IFI35, IFI44, IFI44L
Signaling: IRF1, IRF7, IRF9, STAT1, STAT2
Effectors: RSAD2, BST2, HERC5, USP18, TRIM22

## Antigen Presentation (22 genes)

MHC class I: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, B2M, TAP1, TAP2, TAPBP
MHC class II: HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1
Processing: PSMB8, PSMB9, PSMB10, PSME1, PSME2, CD74, CIITA

## Mito Translation (16 genes)

Elongation: TUFM, GFM1, GFM2, MTIF2
Mitoribosomes: MRPS12, MRPS28, MRPL11, MRPL24, MRPL44
MT aminoacyl-tRNA synthetases: MARS2, KARS1, YARS2, AARS2, EARS2, DARS2, TARS2

Key set for MT-TK (tRNA-Lys) variant 8362T>G — disrupted tRNA charging affects all MT translation.

## mtUPR (13 genes)

MT chaperones: HSPA9, HSPD1, HSPE1, DNAJA3
MT proteases: CLPP, LONP1, AFG3L2, SPG7
Cytosolic stress: HSPB1, HSP90AA1
ISR: ATF4, ATF5, DDIT3

## Glycolysis (16 genes)

Pathway: HK1, HK2, GPI, PFKL, PFKM, ALDOA, TPI1, GAPDH, PGK1, PGAM1, ENO1, PKM, LDHA, LDHB
Transporters: SLC2A1, SLC2A3

Compensatory upregulation expected when OXPHOS is impaired.

## Oxidative Stress (15 genes)

SODs: SOD1, SOD2
Peroxidases: CAT, GPX1, GPX4
Peroxiredoxins/thioredoxins: PRDX1, PRDX2, TXN, TXNRD1
NRF2 pathway: NFE2L2, NQO1, HMOX1
Glutathione: GSR, GCLC, GCLM

## NF-κB Signaling (12 genes)

Subunits: NFKB1, NFKB2, RELA, RELB
IκB/IKK: NFKBIA, IKBKG, CHUK, IKBKB
Targets/feedback: TNFAIP3, BIRC3, TRAF1, TRAF2

## Adding New Gene Sets

1. Define as uppercase character vector: `NEW_GENES <- c("GENE1", ...)`
2. Register in `MODULE_GENE_SETS`:
   ```r
   MODULE_GENE_SETS$new_set <- list(
     genes = NEW_GENES,
     label = "New Set",
     score_col = "new_set_score1"
   )
   ```
3. No other code changes needed — the pipeline automatically iterates over `MODULE_GENE_SETS`
