# Somatic Variants Prediction Analysis Summary

**Date**: 2026-02-17
**Script**: `11.03-variants-annotation-allprediction.R`

## Overview

Analyzed **225 somatic variants** from scMOCHA-data using multiple prediction software tools from two comprehensive databases:
- **MitImpact DB 3.1.3** (APOGEE2 dataset) - 24,190 variants, 152 prediction columns
- **nAPOGEE v1.0.0** - 12,060 variants, focused on non-coding regions

## Prediction Coverage

| Prediction Tool     | Somatic Variants with Predictions | Percentage |
| ------------------- | --------------------------------- | ---------- |
| APOGEE2 (MitImpact) | 107/225                           | 47.6%      |
| nAPOGEE             | 39/225                            | 17.3%      |
| AlphaMissense       | 108/225                           | 48.0%      |
| CADD                | 108/225                           | 48.0%      |

## Prediction Software Categories

### 1. Evolutionary Conservation-based
- **SIFT, SIFT4G**: Sort Intolerant From Tolerant
- **PROVEAN**: Protein Variation Effect Analyzer

### 2. Structural-based
- **PolyPhen2**: Polymorphism Phenotyping v2
- **MutationAssessor**: Functional impact from protein structure

### 3. Machine Learning Ensemble (32 scores in APOGEE2)
- **APOGEE2**: Mitochondrial-specific pathogenicity predictor
- **nAPOGEE**: Non-coding variant predictor
- **AlphaMissense**: DeepMind AI prediction
- **CAROL, Condel**: Combined predictors
- **Meta-SNP**: Meta-predictor combining multiple tools

### 4. Functional/Phenotype-based
- **CADD**: Combined Annotation Dependent Depletion
- **VEST**: Variant Effect Scoring Tool
- **MutationTaster**: Disease-causing potential
- **FATHMM**: Functional Analysis through Hidden Markov Models

### 5. Mitochondria-specific
- **MitoTip**: Mitochondrial variant pathogenicity
- **SNPDryad**: Mitochondrial-specific predictor
- **Mitoclass1**: Mitochondrial variant classifier

## APOGEE2 Classification Distribution

| Classification                                            | Count | Percentage |
| --------------------------------------------------------- | ----- | ---------- |
| Likely-benign                                             | 27    | 25.0%      |
| VUS+ (Variant of Uncertain Significance, lean pathogenic) | 26    | 24.1%      |
| Benign                                                    | 19    | 17.6%      |
| Likely-pathogenic                                         | 12    | 11.1%      |
| VUS- (lean benign)                                        | 12    | 11.1%      |
| VUS                                                       | 10    | 9.3%       |
| Pathogenic                                                | 2     | 1.9%       |

## Top Pathogenic/Likely-Pathogenic Variants

Found **12 unique variants** classified as Pathogenic or Likely-pathogenic by APOGEE2:

| Variant  | Gene    | APOGEE2 Score | APOGEE2 Classification | AlphaMissense     | CADD | PolyPhen2         | SIFT        |
| -------- | ------- | ------------- | ---------------------- | ----------------- | ---- | ----------------- | ----------- |
| 3667T>G  | MT-ND1  | 0.919         | Pathogenic             | likely_pathogenic | 23.5 | probably_damaging | neutral     |
| 9044T>C  | MT-ATP6 | 0.893         | Likely-pathogenic      | likely_pathogenic | 23.6 | probably_damaging | deleterious |
| 14600G>T | MT-ND6  | 0.876         | Likely-pathogenic      | likely_pathogenic | 22.9 | probably_damaging | neutral     |
| 11114T>G | MT-ND4  | 0.852         | Likely-pathogenic      | likely_pathogenic | 23.4 | probably_damaging | deleterious |
| 13592C>T | MT-ND5  | 0.835         | Likely-pathogenic      | likely_pathogenic | 23.8 | probably_damaging | neutral     |
| 10756T>C | MT-ND4L | 0.820         | Likely-pathogenic      | likely_pathogenic | 23.9 | probably_damaging | deleterious |
| 12626C>T | MT-ND5  | 0.794         | Likely-pathogenic      | likely_pathogenic | 23.4 | probably_damaging | neutral     |
| 6967G>A  | MT-CO1  | 0.792         | Likely-pathogenic      | likely_pathogenic | 23.5 | probably_damaging | deleterious |
| 10998G>T | MT-ND4  | 0.784         | Likely-pathogenic      | likely_pathogenic | 23.6 | probably_damaging | deleterious |
| 14487T>G | MT-ND6  | 0.773         | Likely-pathogenic      | likely_pathogenic | 14.3 | probably_damaging | neutral     |
| 13220A>G | MT-ND5  | 0.758         | Likely-pathogenic      | likely_pathogenic | 22.8 | probably_damaging | neutral     |
| 4009A>C  | MT-ND1  | 0.718         | Likely-pathogenic      | likely_pathogenic | 22.7 | probably_damaging | neutral     |

### Gene Distribution of Pathogenic Variants

| Gene    | Classification    | Count |
| ------- | ----------------- | ----- |
| MT-ND5  | Likely-pathogenic | 3     |
| MT-ND4  | Likely-pathogenic | 2     |
| MT-ND6  | Likely-pathogenic | 2     |
| MT-ATP6 | Likely-pathogenic | 1     |
| MT-CO1  | Likely-pathogenic | 1     |
| MT-ND1  | Likely-pathogenic | 1     |
| MT-ND1  | Pathogenic        | 1     |
| MT-ND4L | Likely-pathogenic | 1     |

**Dominant genes**:
- **MT-ND5** (NADH dehydrogenase subunit 5): 3 variants
- **MT-ND4** (NADH dehydrogenase subunit 4): 2 variants
- **MT-ND6** (NADH dehydrogenase subunit 6): 2 variants

## Prediction Score Statistics

For all somatic variants with APOGEE2 predictions:

- **Mean APOGEE2 score**: 0.390
- **Median APOGEE2 score**: 0.364
- **Max APOGEE2 score**: 0.919 (variant 3667T>G in MT-ND1)

## Key Observations

1. **Complex I dominance**: Most pathogenic variants affect NADH dehydrogenase genes (MT-ND1, MT-ND4, MT-ND4L, MT-ND5, MT-ND6), which encode Complex I subunits of the electron transport chain.

2. **Consensus across predictors**: The top pathogenic variants show strong agreement across multiple prediction tools:
   - All 12 pathogenic/likely-pathogenic variants are classified as "likely_pathogenic" by AlphaMissense
   - All have high CADD scores (>14, most >22)
   - All classified as "probably_damaging" by PolyPhen2

3. **SIFT discordance**: Interestingly, SIFT classified 7/12 as "neutral" despite other tools predicting pathogenicity, highlighting the importance of using ensemble predictions.

4. **Coverage gap**: 48% of somatic variants have prediction data, suggesting the need for experimental validation of the remaining 52%.

## Output Files

1. **SOMATIC-VARIANTS-WITH-PREDICTIONS.xlsx**: Complete dataset with all predictions merged
2. **SOMATIC-VARIANTS-PATHOGENIC-TOP.xlsx**: Filtered list of pathogenic/likely-pathogenic variants

## Next Steps

1. Prioritize the 12 pathogenic variants for functional validation
2. Investigate MT-ND5 variants (highest count) for disease associations
3. Examine variants without predictions for potential novel discoveries
4. Correlate pathogenic variants with:
   - Disease phenotypes (AD, COVID-19)
   - Cell type specificity
   - Patient outcomes
