# DIAGRAM: three calling arms, their cutoffs, and their counts

Sample: 7,210 cells, `SAMPLE_LABEL` pending. Every number below is read from
`newplots/compare-with-mgatk/tables/`, generated 2026-09-09. Regenerate the
stage before quoting them.

---

## 1. The three arms end to end

```mermaid
flowchart TB
    RAW["Per-cell strand-resolved allele counts<br/>cell.A/T/C/G.txt.gz<br/>7,210 barcodes"]

    RAW --> CELLA["<b>Cell filter</b><br/>mean MT coverage &gt; 10<br/><b>6,724 cells</b> (486 dropped)"]
    RAW --> CELLB["<b>No cell filter</b><br/><b>7,210 cells</b>"]

    CELLA --> CANDA["Candidate alleles<br/>non-reference, seen on both strands<br/><b>25,580</b>"]
    CELLB --> CANDB["Candidate alleles<br/>non-reference, seen on both strands<br/><b>25,746</b>"]

    CANDA --> CONFA["<b>Confident cell</b><br/>fwd alt &ge; 2 AND rev alt &ge; 2"]
    CANDB --> CONFB["<b>Confident cell</b><br/>fwd alt &ge; 2 AND rev alt &ge; 2<br/>AND fwd+rev alt &ge; 10"]

    CONFA --> S1A["n_cells_conf_detected &ge; 3<br/><b>1,247 variants</b>"]
    CONFB --> S1B["n_cells_conf_detected &ge; 3<br/><b>736 variants</b>"]

    S1A --> GATEA{"<b>VMR / strand gate</b><br/>vmr &gt; 0.01<br/>AND strand r &gt; 0.65"}
    GATEA -->|pass| ARMA["<b>ARM 1<br/>Original mgatk</b><br/><b>230 variants</b>"]
    GATEA -->|"fail strand r only: 628"| DROPA["dropped<br/><b>1,017</b>"]
    GATEA -->|"fail both: 319"| DROPA
    GATEA -->|"fail VMR only: 70"| DROPA

    S1B --> ARMB["<b>ARM 2<br/>scMOCHA variant call</b><br/>no VMR filter<br/>no strand filter<br/><b>736 variants</b>"]

    ARMB --> BL{"<b>Position blacklist</b><br/>mis-alignment + RNA editing"}
    BL -->|"blacklisted: 6"| DROPB["dropped<br/><b>520</b>"]
    BL -->|pass| NC{"<b>Cell-count gate</b><br/>&ge; 10 cells with<br/>AF &ge; 0.05 AND depth &ge; 10"}
    NC -->|"too few cells: 514"| DROPB
    NC -->|pass| ARMC["<b>ARM 3<br/>scMOCHA AF&gt;5%</b><br/><b>216 variants</b>"]

    classDef arm fill:#FFFFFF,stroke:#333,stroke-width:2px
    classDef mg fill:#F7DCC0,stroke:#E18727,stroke-width:2px
    classDef sc fill:#CCE0EE,stroke:#0072B5,stroke-width:2px
    classDef sc5 fill:#B3CBDB,stroke:#044A77,stroke-width:2px
    classDef drop fill:#EFEFEF,stroke:#999,stroke-dasharray:4 3

    class CELLA,CANDA,CONFA,S1A,GATEA mg
    class CELLB,CANDB,CONFB,S1B,BL,NC sc
    class ARMA mg
    class ARMB sc
    class ARMC sc5
    class DROPA,DROPB drop
```

Arm 3 is a strict subset of arm 2: the AF&gt;5% rule is applied downstream in
`src/06.1-collect-variants-new.R`, not during variant calling.

---

## 2. Where the criteria differ

```mermaid
flowchart LR
    subgraph STEP["Pipeline step"]
        direction TB
        A1["1. Cell inclusion"]
        A2["2. Confident cell"]
        A3["3. Variant retention"]
        A4["4. Reliability gate"]
    end

    subgraph MG["Original mgatk"]
        direction TB
        B1["mean MT cov &gt; 10<br/>6,724 / 7,210 cells"]
        B2["fwd &ge; 2 AND rev &ge; 2"]
        B3["n_cells_conf &ge; 3<br/>1,247"]
        B4["vmr &gt; 0.01<br/>AND strand r &gt; 0.65<br/><b>230</b>"]
    end

    subgraph SC["scMOCHA variant call"]
        direction TB
        C1["none<br/>7,210 / 7,210 cells"]
        C2["fwd &ge; 2 AND rev &ge; 2<br/>AND fwd+rev &ge; 10"]
        C3["n_cells_conf &ge; 3<br/><b>736</b>"]
        C4["none"]
    end

    subgraph SC5["scMOCHA AF&gt;5%"]
        direction TB
        D1["none"]
        D2["same as scMOCHA call"]
        D3["same as scMOCHA call"]
        D4["blacklist AND<br/>&ge; 10 cells at AF &ge; 0.05<br/>with depth &ge; 10<br/><b>216</b>"]
    end

    A1 -.-> B1 & C1 & D1
    A2 -.-> B2 & C2 & D2
    A3 -.-> B3 & C3 & D3
    A4 -.-> B4 & C4 & D4

    classDef mg fill:#F7DCC0,stroke:#E18727
    classDef sc fill:#CCE0EE,stroke:#0072B5
    classDef sc5 fill:#B3CBDB,stroke:#044A77
    classDef step fill:#FFFFFF,stroke:#666

    class B1,B2,B3,B4 mg
    class C1,C2,C3,C4 sc
    class D1,D2,D3,D4 sc5
    class A1,A2,A3,A4 step
```

---

## 3. Why the two reported sets differ

Comparing arm 1 (230) against arm 3 (216): 45 shared, 185 mgatk only,
171 scMOCHA AF&gt;5% only.

```mermaid
flowchart LR
    SHARED["<b>45</b><br/>reported by both"]

    ONLY5["<b>171</b><br/>scMOCHA AF&gt;5% only"]
    ONLY5 --> R1["strand r &le; 0.65<br/><b>149</b> (87%)"]
    ONLY5 --> R2["vmr &le; 0.01 AND strand r &le; 0.65<br/><b>17</b>"]
    ONLY5 --> R3["vmr &le; 0.01<br/><b>5</b>"]

    ONLYM["<b>185</b><br/>original mgatk only"]
    ONLYM --> R4["&lt; 10 cells at AF &ge; 0.05<br/><b>117</b> (63%)"]
    ONLYM --> R5["&lt; 3 cells with fwd+rev alt &ge; 10<br/><b>68</b> (37%)"]

    classDef mg fill:#F7DCC0,stroke:#E18727
    classDef sc5 fill:#B3CBDB,stroke:#044A77
    classDef sh fill:#DDDDDD,stroke:#666
    classDef reason fill:#FFFFFF,stroke:#999

    class ONLYM mg
    class ONLY5 sc5
    class SHARED sh
    class R1,R2,R3,R4,R5 reason
```

Read this as the summary of the whole stage:

- **87% of what mgatk misses is lost to the strand-correlation cutoff alone**,
  not to VMR and not to read depth.
- **What scMOCHA drops is dropped for prevalence or read support**: too few
  carrier cells, or too few alt reads per cell. The variants only mgatk reports
  carry a median of 9 alt reads per carrier cell, against 48 for the shared
  variants.

---

## 4. Notes

- Every count comes from `03-funnel-counts.tsv`, `03-exclusion-reasons.tsv`,
  `02-cell-inclusion.tsv` and `05-read-support.tsv`.
- Exclusion reasons are assigned by first failed criterion in each arm's own
  order, so each excluded variant is attributed to exactly one cutoff. The
  order is set in `fn_exclusion_*` in `config.R`.
- The mgatk S1 -> S2 breakdown (628 / 319 / 70) counts all 1,247 mgatk S1
  variants, not only those scMOCHA also reports.
- Colours match `color_arm` in `high-res/00-colors.R`: mgatk `#E18727`,
  scMOCHA call `#0072B5`, scMOCHA AF&gt;5% `#044A77`. The fills above are
  lightened versions for readability behind black text.
