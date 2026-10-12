# GI-melanoma-classifier

Code and processed data for:

**Tao H, Zhang J, Musee HK. Integrated Multi-Feature Genomic Instability Profiling Identifies Distinct Melanoma Tumor Phenotypes Associated with Copy Number Alterations. Cancer Informatics (submitted, 2026).**

## Contents

### Analysis scripts

| File | Description |
|---|---|
| `classify_new_samples.R` | Standalone R script to classify new melanoma samples as GI-High or GI-Low using TCGA-derived parameters (only base R required) |
| `Cluster_stability.R` | Bootstrap stability and consensus clustering of the TCGA discovery cohort |
| `Independent_validation.R` | Computation of independent validation features (focal amplifications, focal deletions, genome-wide gain/loss fractions, mutational signature proxy) |

### TCGA discovery cohort results (`TCGA/`)

| File | Description |
|---|---|
| `TCGA_FINAL_TWO_PHENOTYPE_MODEL.csv` | TCGA samples (n = 425) with GI-High / GI-Low assignments |
| `TCGA_WITH_CLUSTERS.csv` | TCGA samples with raw cluster labels |
| `TCGA_CLUSTER_CENTERS.csv` | k-means centroids in standardized feature space |
| `TCGA_Bootstrap_Stability_Results.csv` | Bootstrap ARI metrics (mean ARI = 0.959) |
| `TCGA_Consensus_Clustering_Results.csv` | Consensus clustering PAC score (0.033) |
| `TCGA_TWO_PHENOTYPE_SUMMARY.csv` | Phenotype summary statistics |
| `Cluster_Summary_Final.csv` | Cluster centroids |
| `Gap_Statistic.csv` | Gap statistic output |
| `Cox_Two_Phenotype_Model.txt` | Cox survival model output |
| `Cox_Final_Model.txt` | Full Cox model output |
| `*.png` | Analysis figures |

### MSK-IMPACT validation cohort results (`MSK/`)

| File | Description |
|---|---|
| `MSK_VALIDATED_PHENOTYPES.csv` | MSK-IMPACT samples (n = 690) with GI-High / GI-Low assignments |
| `MSK_Independent_Validation_Features.csv` | Independent validation features (focal events, gain/loss fractions, C>T fraction) |
| `MSK_Independent_Validation_Summary.csv` | Summary statistics for independent validation |
| `MSK_Chromosome_CNA_Frequency_FINAL.csv` | Chromosome-level CNA frequencies (Table 4 source) |
| `MSK_Chromosome_CNA_Fisher_FINAL.csv` | Fisher's exact test results for chromosome CNA |
| `MSK_Driver_Mutation_Frequency.csv` | Driver mutation frequencies (Table 3 source) |
| `MSK_Driver_Mutation_Fisher_Test.csv` | Fisher's exact test results for driver mutations |
| `MSK_PHENOTYPE_SUMMARY.csv` | MSK phenotype summary statistics |
| `MSK_Effect_Size.csv` | Fold change calculations |
| `MSK_Wilcoxon_Feature_Comparison.csv` | Wilcoxon test results |
| `MSK_Clinical_Overview_By_Phenotype.csv` | Clinical summary by phenotype |
| `MSK_Clinical_Continuous_Comparison.csv` | Continuous clinical comparisons |
| `MSK_Clinical_Categorical_Comparison.csv` | Categorical clinical comparisons |
| `MSK_FGA_Aneuploidy_Correlation.txt` | Spearman correlation between FGA and aneuploidy |
| `*.png` | Analysis figures |

## How to classify new samples

```r
source("classify_new_samples.R")

# Prepare a data frame with three columns: FGA, Aneuploidy_Arm_Score, logTMB
new_samples <- data.frame(
  FGA                  = c(0.05, 0.45, 0.62),
  Aneuploidy_Arm_Score = c(2,    18,   27),
  logTMB               = c(3.2,  6.1,  6.8)
)

result <- classify_GI(new_samples)
print(result)
