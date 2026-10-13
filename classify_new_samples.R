library(data.table)

TCGA_RESULT <- "C:/Users/Henry/Desktop/SeerStat/Seer_Rstat/Results_TCGA"

# Load TCGA centroids and data
centers <- fread(file.path(TCGA_RESULT, "TCGA_CLUSTER_CENTERS.csv"))
tcga    <- fread(file.path(TCGA_RESULT, "TCGA_FINAL_TWO_PHENOTYPE_MODEL.csv"))

# TCGA-derived means and SDs (used for MSK scaling)
means <- apply(tcga[, c("FGA", "Aneuploidy_Arm_Score", "log_mutation_count")], 2, mean)
sds   <- apply(tcga[, c("FGA", "Aneuploidy_Arm_Score", "log_mutation_count")], 2, sd)

cat("TCGA-DERIVED SCALING PARAMETERS")
cat("\nMeans (center):\n")
print(round(means, 6))
cat("\nStandard deviations (scale):\n")
print(round(sds, 6))

cat("TCGA K-MEANS CLUSTER CENTROIDS (in z-score space)")
print(round(centers, 6))

# Also identify which centroid is GI-High
centroid_high <- which.max(centers$FGA)
cat("\nGI-High is centroid row:", centroid_high, "\n")
cat("GI-Low  is centroid row:", 3 - centroid_high, "\n")





#
# classify_new_samples.R
#
# Purpose
#
# Assign new melanoma samples to one of two genomic instability phenotypes
# (GI-High or GI-Low) using the k-means classifier that was trained on the
# TCGA Skin Cutaneous Melanoma (SKCM) discovery cohort (n = 425).
#
# This script reproduces the external validation procedure used in:
#
#   Tao H, Zhang J, Musee HK. Integrated Multi-Feature Genomic Instability
#   Profiling Identifies Distinct Melanoma Tumor Phenotypes Associated with
#   Copy Number Alterations. Cancer Informatics (submitted, 2026).
#
# It is intentionally self-contained: it depends only on base R, so it can be
# executed by any researcher without installing additional packages.
#
# Input
# -----
# A data.frame (or matrix coercible to a data.frame) with at least three
# numeric columns, named EXACTLY:
#
#   FGA                   Fraction of genome altered (0 to 1)
#   Aneuploidy_Arm_Score  Number of chromosome arms altered (integer)
#   log_mutation_count    log(1 + nonsynonymous mutation count)
#
# Output
# 
# The original data.frame with four additional columns:
#
#   Cluster         Integer cluster label (1 = GI-High, 2 = GI-Low)
#   Phenotype       Factor with levels "GI-Low" and "GI-High"
#   dist_GI_High    Squared Euclidean distance to the GI-High centroid
#   dist_GI_Low     Squared Euclidean distance to the GI-Low centroid
#
# Usage
# 
#   source("classify_new_samples.R")
#
#   my_samples <- data.frame(
#     FGA                  = c(0.05, 0.45, 0.62),
#     Aneuploidy_Arm_Score = c(2,    18,   27),
#     log_mutation_count   = c(3.2,  6.1,  6.8)
#   )
#
#   classify_GI(my_samples)
#
# 
# SECTION 1. TCGA-DERIVED PARAMETERS (DO NOT MODIFY)
#
# These numbers are the means, standard deviations, and k-means centroids
# extracted from the TCGA discovery cohort. They must remain unchanged so
# that new samples are classified on exactly the same scale as the original
# MSK validation.

# Feature means in the TCGA discovery cohort (used for centering)
TCGA_MEANS <- c(
  FGA                  = 0.319367,
  Aneuploidy_Arm_Score = 11.117647,
  log_mutation_count   = 5.983047
)

# Feature standard deviations in the TCGA discovery cohort (used for scaling)
TCGA_SDS <- c(
  FGA                  = 0.206321,
  Aneuploidy_Arm_Score = 8.220734,
  log_mutation_count   = 1.192553
)

# k-means cluster centroids in standardized (z-score) space.
#
# The two rows correspond to the two clusters identified by k-means on the
# TCGA discovery cohort. Row 1 is the cluster whose centroid has the higher
# FGA value in standardized space and is therefore labeled GI-High. Row 2 is
# labeled GI-Low.
TCGA_CENTROIDS <- matrix(
  data = c(
    #  FGA        Aneuploidy_Arm_Score   log_mutation_count
    1.250271,       1.278153,        0.066189,   # Row 1 = GI-High
    -0.538837,      -0.550854,       -0.028526    # Row 2 = GI-Low
  ),
  nrow     = 2,
  byrow    = TRUE,
  dimnames = list(
    c("GI-High", "GI-Low"),
    c("FGA", "Aneuploidy_Arm_Score", "log_mutation_count")
  )
)

# The order of the columns used throughout the classifier. This ensures that
# the mean vector, SD vector, and centroid matrix are always aligned.
FEATURE_ORDER <- c("FGA", "Aneuploidy_Arm_Score", "log_mutation_count")


# 
# SECTION 2. THE CLASSIFIER
# 
#' Classify melanoma samples as GI-High or GI-Low
#'
#' @param df A data.frame with columns FGA, Aneuploidy_Arm_Score, logTMB.
#' @return The input data.frame with Cluster, Phenotype, dist_GI_High, and
#'         dist_GI_Low appended.
#' @export
classify_GI <- function(df) {
  
  # 2.1 Validate input 
  
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  
  missing_cols <- setdiff(FEATURE_ORDER, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Input data.frame is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ". Required columns are: ",
      paste(FEATURE_ORDER, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  #  2.2 Extract and validate feature matrix
  
  X <- as.matrix(df[, FEATURE_ORDER, drop = FALSE])
  storage.mode(X) <- "double"
  
  if (anyNA(X)) {
    stop(
      "Input contains missing values in one or more feature columns. ",
      "Please remove or impute missing values before classification.",
      call. = FALSE
    )
  }
  
  #2.3 Standardize using TCGA-derived means and SDs 
  #
  # For each feature j:
  #     z_j = (x_j - mean_j) / sd_j
  #
  # Note: we deliberately use the TCGA-derived mean and SD, not the mean and
  # SD of the new samples, so that the new samples are placed on the same
  # scale as the original TCGA training data.
  
  X_centered <- sweep(X, MARGIN = 2, STATS = TCGA_MEANS[FEATURE_ORDER], FUN = "-")
  X_scaled   <- sweep(X_centered, MARGIN = 2, STATS = TCGA_SDS[FEATURE_ORDER], FUN = "/")
  
  # 2.4 Compute squared Euclidean distance to each centroid 
  #
  # For each sample i and centroid k:
  #     d^2(i, k) = sum_j ( z_ij - c_kj )^2
  
  n_samples <- nrow(X_scaled)
  n_features <- length(FEATURE_ORDER)
  
  # Replicate the centroid rows to match the number of samples for vectorized
  # computation. Each centroid becomes an n_samples x n_features matrix.
  high_matrix <- matrix(
    TCGA_CENTROIDS["GI-High", ],
    nrow = n_samples, ncol = n_features, byrow = TRUE
  )
  low_matrix <- matrix(
    TCGA_CENTROIDS["GI-Low", ],
    nrow = n_samples, ncol = n_features, byrow = TRUE
  )
  
  dist_high <- rowSums((X_scaled - high_matrix)^2)
  dist_low  <- rowSums((X_scaled - low_matrix)^2)
  
  # 2.5 Assign each sample to the nearest centroid 
  
  phenotype <- ifelse(dist_high < dist_low, "GI-High", "GI-Low")
  cluster   <- ifelse(phenotype == "GI-High", 1L, 2L)
  
  #  2.6 Return augmented data.frame 
  out <- df
  out$Cluster      <- cluster
  out$Phenotype    <- factor(phenotype, levels = c("GI-Low", "GI-High"))
  out$dist_GI_High <- dist_high
  out$dist_GI_Low  <- dist_low
  out
}


# # SECTION 3. WORKED EXAMPLE (for documentation purposes)
#
# Running this block reproduces a small demonstration. It is wrapped in
# if (FALSE) so it does not execute when the file is sourced.
# 

if (FALSE) {
  
  # Three hypothetical samples
  demo <- data.frame(
    FGA                  = c(0.05, 0.45, 0.62),
    Aneuploidy_Arm_Score = c(2,    18,   27),
    log_mutation_count   = c(3.2,  6.1,  6.8)
  )
  
  # Classify
  result <- classify_GI(demo)
  
  # Inspect output
  print(result[, c("FGA", "Aneuploidy_Arm_Score", "log_mutation_count",
                   "Phenotype", "dist_GI_High", "dist_GI_Low")])
  
  # Expected: sample 1 (low FGA, low aneuploidy) -> GI-Low
  #           sample 2 (mid FGA, mid aneuploidy) -> intermediate
  #           sample 3 (high FGA, high aneuploidy) -> GI-High
}


# 
# SECTION 4. REPRODUCIBILITY NOTE
#
# The classifier above reproduces the exact external validation procedure
# described in the Methods section of the manuscript. Given the MSK-IMPACT
# feature table (columns FGA, Aneuploidy_Arm_Score, logTMB), running
# classify_GI() returns the same GI-High / GI-Low assignments that were used
# in the manuscript.
# 
