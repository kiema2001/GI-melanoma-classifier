# 
#figures.R
#
# 

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(survival)
library(survminer)
library(ggpubr)

TCGA_RESULT <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/Results_TCGA"
MSK_RESULT  <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/Results_MSK"
FIG_DIR     <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/Figures_FINAL"
dir.create(FIG_DIR, showWarnings = FALSE)

phenotype_colors <- c("GI-Low" = "#4DBBD5", "GI-High" = "#E64B35")

tcga <- fread(file.path(TCGA_RESULT, "TCGA_FINAL_TWO_PHENOTYPE_MODEL.csv"),
              data.table = FALSE)
msk  <- fread(file.path(MSK_RESULT, "MSK_VALIDATED_PHENOTYPES.csv"),
              data.table = FALSE)
msk$Phenotype <- factor(msk$Phenotype, levels = c("GI-Low", "GI-High"))

# Figure 1
pca <- prcomp(scale(tcga[, c("FGA", "Aneuploidy_Arm_Score",
                             "log_mutation_count")]))
variance <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Phenotype = tcga$Phenotype)

p1A <- ggplot(pca_df, aes(PC1, PC2, color = Phenotype)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = phenotype_colors) +
  labs(title = "A. Principal component analysis",
       x = paste0("PC1 (", variance[1], "%)"),
       y = paste0("PC2 (", variance[2], "%)")) +
  theme_classic(base_size = 13)

sil_df <- data.frame(k = 2:6,
                     Silhouette = c(0.320, 0.316, 0.303, 0.290, 0.297))

p1B <- ggplot(sil_df, aes(k, Silhouette)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  labs(title = "B. Cluster optimization",
       x = "Number of clusters", y = "Average silhouette width") +
  theme_classic(base_size = 13)

p1C <- ggplot(tcga, aes(Phenotype, FGA, fill = Phenotype)) +
  geom_boxplot() + scale_fill_manual(values = phenotype_colors) +
  labs(title = "C. Fraction of genome altered", y = "FGA", x = "Phenotype") +
  theme_classic(base_size = 13) + theme(legend.position = "none")

figure1 <- (p1A + p1B) / p1C
ggsave(file.path(FIG_DIR, "Figure1_TCGA_Discovery.png"),
       figure1, width = 10, height = 9, dpi = 300)

# Figure 2
feature_long <- tcga %>%
  select(Phenotype, FGA, Aneuploidy_Arm_Score, log_mutation_count) %>%
  pivot_longer(cols = -Phenotype, names_to = "Feature", values_to = "Value")

figure2 <- ggplot(feature_long, aes(Phenotype, Value, fill = Phenotype)) +
  geom_boxplot() + facet_wrap(~Feature, scales = "free") +
  scale_fill_manual(values = phenotype_colors) +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "TCGA genomic characteristics") +
  theme_classic(base_size = 13) + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "Figure2_TCGA_Genomic_Features.png"),
       figure2, width = 10, height = 5, dpi = 300)

# Figure 3
km_fit <- survfit(Surv(OS_time, OS_event) ~ Phenotype, data = tcga)
km_plot <- ggsurvplot(km_fit, data = tcga,
                      risk.table = TRUE, conf.int = TRUE, pval = TRUE,
                      palette = phenotype_colors,
                      legend.title = "Genomic instability phenotype",
                      xlab = "Overall survival (months)",
                      ylab = "Survival probability",
                      ggtheme = theme_classic())
ggsave(file.path(FIG_DIR, "Figure3_TCGA_Survival.png"),
       km_plot$plot, width = 8, height = 6, dpi = 300)

# Figure 4 — MSK validation
msk_pca <- prcomp(scale(msk[, c("FGA", "Aneuploidy_Arm_Score",
                                "log_mutation_count")]))
msk_var <- round(100 * (msk_pca$sdev^2 / sum(msk_pca$sdev^2))[1:2], 1)
pca_msk <- data.frame(PC1 = msk_pca$x[, 1], PC2 = msk_pca$x[, 2],
                      Phenotype = msk$Phenotype)

p4A <- ggplot(pca_msk, aes(PC1, PC2, color = Phenotype)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = phenotype_colors) +
  labs(title = "A. MSK PCA (confirmatory)",
       x = paste0("PC1 (", msk_var[1], "%)"),
       y = paste0("PC2 (", msk_var[2], "%)")) +
  theme_classic(base_size = 12) + theme(legend.position = "right")

msk_indep <- fread(file.path(MSK_RESULT,
                             "MSK_Independent_Validation_Features.csv"),
                   data.table = FALSE)
msk_indep$Phenotype <- factor(msk_indep$Phenotype,
                              levels = c("GI-Low", "GI-High"))

make_box <- function(df, feature, title) {
  ggplot(df, aes(Phenotype, .data[[feature]], fill = Phenotype)) +
    geom_boxplot(outlier.size = 0.4, width = 0.6) +
    scale_fill_manual(values = phenotype_colors) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       size = 3.2, label.y.npc = 0.95) +
    labs(title = title, y = "Value") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.title.x = element_blank())
}

p4B <- make_box(msk_indep, "n_focal_amp", "B. Focal amplifications")
p4C <- make_box(msk_indep, "n_focal_del", "C. Focal deletions")
p4D <- make_box(msk_indep, "frac_loss",   "D. Genome loss fraction")
p4E <- make_box(msk_indep, "frac_gain",   "E. Genome gain fraction")
p4F <- make_box(msk_indep, "frac_CT",     "F. C>T fraction (control)")

figure4 <- (p4A | p4B | p4C) / (p4D | p4E | p4F)
ggsave(file.path(FIG_DIR, "Figure4_MSK_External_Validation.png"),
       figure4, width = 14, height = 9, dpi = 300)

# Figure 5
mut_freq <- fread(file.path(MSK_RESULT, "MSK_Driver_Mutation_Frequency.csv"))
fisher   <- fread(file.path(MSK_RESULT, "MSK_Driver_Mutation_Fisher.csv"))

p5A <- ggplot(mut_freq, aes(Hugo_Symbol, Frequency_percent, fill = Phenotype)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = phenotype_colors) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "A. Driver mutation frequency",
       x = "Driver gene", y = "Mutated samples (%)")

p5B <- ggplot(fisher, aes(Odds_Ratio, reorder(Gene, Odds_Ratio))) +
  geom_point(size = 3) + geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "B. Odds ratio", x = "Odds ratio", y = "Gene") +
  theme_classic(base_size = 13)

p5C <- ggplot(fisher, aes(Gene, -log10(FDR))) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "C. FDR significance", x = "Gene", y = "-log10(FDR)")

figure5 <- p5A / (p5B | p5C)
ggsave(file.path(FIG_DIR, "Figure5_Driver_Mutations.png"),
       figure5, width = 12, height = 10, dpi = 300)

# Figure 6
cna_freq <- fread(file.path(MSK_RESULT,
                            "MSK_Chromosome_CNA_Frequency_FINAL.csv"))
cna_freq <- cna_freq %>% filter(!is.na(Phenotype))

p6A <- ggplot(cna_freq, aes(factor(Chromosome), Frequency, fill = Phenotype)) +
  geom_col(position = "dodge") +
  facet_wrap(~dominant_call) +
  scale_fill_manual(values = phenotype_colors) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  labs(title = "A. Chromosome-level CNA frequency",
       x = "Chromosome", y = "Altered samples (%)")

cna_table <- cna_freq %>%
  select(Chromosome, dominant_call, Phenotype, Patients) %>%
  pivot_wider(names_from = Phenotype, values_from = Patients, values_fill = 0) %>%
  mutate(Odds_Ratio = (`GI-High` / sum(`GI-High`)) /
           (`GI-Low`  / sum(`GI-Low`)))

p6B <- ggplot(cna_table, aes(factor(Chromosome), Odds_Ratio,
                             color = dominant_call)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~dominant_call) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "B. CNA enrichment in GI-High tumors",
       x = "Chromosome", y = "Relative enrichment")

figure6 <- p6A / p6B
ggsave(file.path(FIG_DIR, "Figure6_Chromosome_CNA.png"),
       figure6, width = 12, height = 10, dpi = 300)

