# 
# 01_TCGA_discovery.R
#
# Purpose
# 
# Reproduce the TCGA SKCM discovery analysis:
#   - Compute FGA, arm-level aneuploidy score, and log-transformed
#     nonsynonymous mutation count
#   - k-means clustering (k = 2)
#   - Cox proportional hazards regression
#   - Save processed outputs
#
# Inputs
# 
#   TCGA_DIR/data_mutations.txt
#   TCGA_DIR/data_cna_hg19.seg
#   TCGA_DIR/data_clinical_patient.txt
#   TCGA_DIR/data_clinical_sample.txt
#
# Outputs (in TCGA_RESULT)
# 
#   TCGA_FINAL_TWO_PHENOTYPE_MODEL.csv
#   TCGA_CLUSTER_CENTERS.csv
#   TCGA_TWO_PHENOTYPE_SUMMARY.csv
#   Cox_Two_Phenotype_Model.txt

#load librariess

library(data.table)
library(dplyr)
library(tidyr)
library(survival)

# Paths (modify to your local setup)
TCGA_DIR    <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/skcm_tcga_pan_can_atlas_2018"
TCGA_RESULT <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/Results_TCGA"
dir.create(TCGA_RESULT, showWarnings = FALSE)

# 1. Load raw data 
clinical_patient <- fread(file.path(TCGA_DIR, "data_clinical_patient.txt"),
                          skip = 4, data.table = FALSE)
clinical_sample  <- fread(file.path(TCGA_DIR, "data_clinical_sample.txt"),
                          skip = 4, data.table = FALSE)
seg              <- fread(file.path(TCGA_DIR, "data_cna_hg19.seg"),
                          data.table = FALSE)
mut              <- fread(file.path(TCGA_DIR, "data_mutations.txt"),
                          data.table = FALSE)

# 2. Compute FGA per sample 
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

# 3. Compute arm-level aneuploidy score 
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

#  4. Compute nonsynonymous mutation count 
nonsyn_classes <- c("Missense_Mutation", "Nonsense_Mutation",
                    "Frame_Shift_Del", "Frame_Shift_Ins",
                    "Splice_Site", "Translation_Start_Site")

TMB <- mut %>%
  filter(Variant_Classification %in% nonsyn_classes) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(NonSyn_Mutations = n(), .groups = "drop") %>%
  mutate(logTMB = log1p(NonSyn_Mutations))

#5. Merge 
tcga <- clinical_sample %>%
  left_join(clinical_patient, by = "PATIENT_ID") %>%
  left_join(FGA, by = c("SAMPLE_ID" = "ID")) %>%
  left_join(aneuploidy, by = c("SAMPLE_ID" = "ID")) %>%
  left_join(TMB, by = c("SAMPLE_ID" = "Tumor_Sample_Barcode")) %>%
  drop_na(FGA, Aneuploidy_Arm_Score, logTMB) %>%
  mutate(
    OS_time  = as.numeric(OS_MONTHS),
    OS_event = ifelse(grepl("DECEASED", OS_STATUS), 1, 0),
    AGE      = as.numeric(AGE),
    SEX      = factor(SEX)
  ) %>%
  drop_na(OS_time, OS_event)

cat("TCGA cohort n =", nrow(tcga), "\n")

#  6. k-means clustering (k = 2)
features <- tcga %>% select(FGA, Aneuploidy_Arm_Score, logTMB)
scaled   <- scale(features)

set.seed(123)
km <- kmeans(scaled, centers = 2, nstart = 100, iter.max = 100)

# Label GI-High as the cluster with the higher FGA centroid
centers <- as.data.frame(km$centers)
high_cluster <- which.max(centers$FGA)

tcga$Cluster <- factor(km$cluster)
tcga$Phenotype <- ifelse(tcga$Cluster == high_cluster, "GI-High", "GI-Low")
tcga$Phenotype <- factor(tcga$Phenotype, levels = c("GI-Low", "GI-High"))

cat("Phenotype distribution:\n")
print(table(tcga$Phenotype))

# 7. Cox regression
cox <- coxph(Surv(OS_time, OS_event) ~ Phenotype + AGE + SEX + logTMB,
             data = tcga)
print(summary(cox))

# 8. Save outputs
write.csv(tcga,
          file.path(TCGA_RESULT, "TCGA_FINAL_TWO_PHENOTYPE_MODEL.csv"),
          row.names = FALSE)
write.csv(as.data.frame(km$centers),
          file.path(TCGA_RESULT, "TCGA_CLUSTER_CENTERS.csv"),
          row.names = FALSE)
capture.output(summary(cox),
               file = file.path(TCGA_RESULT, "Cox_Two_Phenotype_Model.txt"))

summary_table <- tcga %>%
  group_by(Phenotype) %>%
  summarise(n = n(),
            FGA = mean(FGA),
            Aneuploidy = mean(Aneuploidy_Arm_Score),
            logTMB = mean(logTMB),
            .groups = "drop")
write.csv(summary_table,
          file.path(TCGA_RESULT, "TCGA_TWO_PHENOTYPE_SUMMARY.csv"),
          row.names = FALSE)

cat("\nTCGA discovery analysis complete.\n")
