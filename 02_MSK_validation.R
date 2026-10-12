# 02_MSK_validation.R
#
# Purpose
#
# Reproduce the MSK-IMPACT external validation:
#   - Compute FGA, arm-level aneuploidy score, and log-transformed
#     nonsynonymous mutation count
#   - Classify samples using TCGA-derived centroids
#   - Compare classifier features between phenotypes (reproducibility)
#
# Inputs
#
#   MSK_DIR/data_mutations.txt
#   MSK_DIR/data_cna_hg19.seg
#   MSK_DIR/data_clinical_patient.txt
#   MSK_DIR/data_clinical_sample.txt
#   TCGA_RESULT/TCGA_CLUSTER_CENTERS.csv
#   TCGA_RESULT/TCGA_FINAL_TWO_PHENOTYPE_MODEL.csv
#
# Outputs (in MSK_RESULT)
# 
#   MSK_VALIDATED_PHENOTYPES.csv
#   MSK_PHENOTYPE_SUMMARY.csv
#   MSK_Wilcoxon_Feature_Comparison.csv
#   MSK_Effect_Size.csv
#

library(data.table)
library(dplyr)
library(tidyr)

MSK_DIR     <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/mel_mskimpact_2020"
TCGA_RESULT <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/Results_TCGA"
MSK_RESULT  <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/Results_MSK"
dir.create(MSK_RESULT, showWarnings = FALSE)

#  1. Load raw data 
patient <- fread(file.path(MSK_DIR, "data_clinical_patient.txt"),
                 skip = 4, data.table = FALSE)
sample  <- fread(file.path(MSK_DIR, "data_clinical_sample.txt"),
                 skip = 4, data.table = FALSE)
seg     <- fread(file.path(MSK_DIR, "data_cna_hg19.seg"),
                 data.table = FALSE)
mut     <- fread(file.path(MSK_DIR, "data_mutations.txt"),
                 data.table = FALSE)

#2. FGA 
seg$segment_length <- seg$loc.end - seg$loc.start + 1
threshold <- 0.20

FGA <- seg %>%
  mutate(altered = abs(seg.mean) >= threshold) %>%
  group_by(ID) %>%
  summarise(
    altered_bp = sum(segment_length[altered], na.rm = TRUE),
    total_bp   = sum(segment_length, na.rm = TRUE),
    FGA        = altered_bp / total_bp,
    .groups    = "drop"
  )

#  3. Arm-level aneuploidy 
centromeres <- data.frame(
  chrom = 1:22,
  cent  = c(125000000, 93300000, 91000000, 50000000, 48800000,
            61000000, 59900000, 45000000, 49000000, 40200000,
            53700000, 35800000, 17900000, 17600000, 19000000,
            35300000, 22200000, 15400000, 24600000, 26300000,
            12000000, 15000000)
)

seg_clean <- seg %>%
  mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>%
  filter(chrom %in% 1:22) %>%
  left_join(centromeres, by = "chrom") %>%
  mutate(arm = ifelse(loc.end <= cent,
                      paste0(chrom, "p"),
                      paste0(chrom, "q")))

arm_scores <- seg_clean %>%
  mutate(altered = abs(seg.mean) >= threshold) %>%
  group_by(ID, arm) %>%
  summarise(
    altered_fraction = sum(segment_length[altered], na.rm = TRUE) /
      sum(segment_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(arm_altered = altered_fraction >= 0.8)

aneuploidy <- arm_scores %>%
  group_by(ID) %>%
  summarise(Aneuploidy_Arm_Score = sum(arm_altered, na.rm = TRUE),
            .groups = "drop")

# 4. Nonsynonymous mutation count 
nonsyn_classes <- c("Missense_Mutation", "Nonsense_Mutation",
                    "Frame_Shift_Del", "Frame_Shift_Ins",
                    "Splice_Site", "Translation_Start_Site")

TMB <- mut %>%
  filter(Variant_Classification %in% nonsyn_classes) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(NonSyn_Mutations = n(), .groups = "drop") %>%
  mutate(logTMB = log1p(NonSyn_Mutations))

#  5. Merge MSK data 
msk <- sample %>%
  left_join(patient, by = "PATIENT_ID") %>%
  left_join(FGA, by = c("SAMPLE_ID" = "ID")) %>%
  left_join(aneuploidy, by = c("SAMPLE_ID" = "ID")) %>%
  left_join(TMB, by = c("SAMPLE_ID" = "Tumor_Sample_Barcode")) %>%
  mutate(
    NonSyn_Mutations = ifelse(is.na(NonSyn_Mutations), 0, NonSyn_Mutations),
    logTMB           = ifelse(is.na(logTMB), 0, logTMB)
  ) %>%
  drop_na(FGA, Aneuploidy_Arm_Score, logTMB)

cat("MSK cohort n =", nrow(msk), "\n")

#6. Classify MSK using TCGA-derived centroids
centers <- fread(file.path(TCGA_RESULT, "TCGA_CLUSTER_CENTERS.csv"))
tcga    <- fread(file.path(TCGA_RESULT, "TCGA_FINAL_TWO_PHENOTYPE_MODEL.csv"))

means <- apply(tcga[, c("FGA", "Aneuploidy_Arm_Score", "logTMB")], 2, mean)
sds   <- apply(tcga[, c("FGA", "Aneuploidy_Arm_Score", "logTMB")], 2, sd)

msk_scaled <- scale(msk[, c("FGA", "Aneuploidy_Arm_Score", "logTMB")],
                    center = means, scale = sds)

dist_to <- function(x, centroid) sum((x - centroid)^2)

clusters <- apply(msk_scaled, 1, function(x) {
  d1 <- dist_to(x, as.numeric(centers[1, ]))
  d2 <- dist_to(x, as.numeric(centers[2, ]))
  ifelse(d1 < d2, 1, 2)
})

msk$Cluster <- factor(clusters)
centroid_high <- which.max(centers$FGA)
msk$Phenotype <- ifelse(msk$Cluster == centroid_high, "GI-High", "GI-Low")
msk$Phenotype <- factor(msk$Phenotype, levels = c("GI-Low", "GI-High"))

cat("MSK phenotype distribution:\n")
print(table(msk$Phenotype))

# 7. Save 
write.csv(msk,
          file.path(MSK_RESULT, "MSK_VALIDATED_PHENOTYPES.csv"),
          row.names = FALSE)

#8. Summary 
summary_table <- msk %>%
  group_by(Phenotype) %>%
  summarise(n = n(),
            FGA = mean(FGA),
            Aneuploidy = mean(Aneuploidy_Arm_Score),
            logTMB = mean(logTMB),
            .groups = "drop")
print(summary_table)
write.csv(summary_table,
          file.path(MSK_RESULT, "MSK_PHENOTYPE_SUMMARY.csv"),
          row.names = FALSE)

#9. Wilcoxon tests
wilcox_results <- data.frame()
for (v in c("FGA", "Aneuploidy_Arm_Score", "logTMB")) {
  test <- wilcox.test(msk[[v]] ~ msk$Phenotype)
  wilcox_results <- rbind(wilcox_results, data.frame(
    Feature = v, Statistic = test$statistic, P_value = test$p.value
  ))
}
print(wilcox_results)
write.csv(wilcox_results,
          file.path(MSK_RESULT, "MSK_Wilcoxon_Feature_Comparison.csv"),
          row.names = FALSE)

# 10. Effect sizes
effect_results <- data.frame()
for (v in c("FGA", "Aneuploidy_Arm_Score", "logTMB")) {
  x <- msk %>% filter(Phenotype == "GI-High") %>% pull(all_of(v))
  y <- msk %>% filter(Phenotype == "GI-Low")  %>% pull(all_of(v))
  effect_results <- rbind(effect_results, data.frame(
    Feature = v,
    GI_High_Mean = mean(x),
    GI_Low_Mean  = mean(y),
    Fold_Change  = mean(x) / mean(y)
  ))
}
print(effect_results)
write.csv(effect_results,
          file.path(MSK_RESULT, "MSK_Effect_Size.csv"),
          row.names = FALSE)

cat("\nMSK validation complete.\n")