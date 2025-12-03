#' -----------------------------------------------------------------------------
#' Data Preparation Script: Combined Human and Mouse TF-Target Interactions
#' -----------------------------------------------------------------------------
#'
#' Source: DoRothEA (Discriminant Regulon Expression Analysis)
#' Purpose: To retrieve, clean, and combine high-confidence transcription factor (TF)
#'          to target gene interactions for Human and Mouse into a single dataset
#'          for the TCMDATA package.
#'
#' References:
#' 1. https://saezlab.github.io/dorothea/
#' 2. Garcia-Alonso L, et al. Genome Research. 2019. DOI: 10.1101/gr.240663.118.
#' -----------------------------------------------------------------------------

# 1. Load necessary libraries
library(yulab.utils)
pload(dorothea)
pload(dplyr)
pload(usethis)

# -----------------------------------------------------------------------------
# Part A: Process Human Data
# -----------------------------------------------------------------------------
data(dorothea_hs, package = "dorothea")

human_data <- dorothea_hs %>%
  # Filter for high-to-medium confidence interactions (A, B, C)
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  # Select necessary columns
  dplyr::select(tf, target, confidence, mor) %>%
  # Add Species column
  dplyr::mutate(Species = "Human")

# -----------------------------------------------------------------------------
# Part B: Process Mouse Data
# -----------------------------------------------------------------------------
data(dorothea_mm, package = "dorothea")

mouse_data <- dorothea_mm %>%
  # Filter for high-to-medium confidence interactions (A, B, C)
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  # Select necessary columns
  dplyr::select(tf, target, confidence, mor) %>%
  # Add Species column
  dplyr::mutate(Species = "Mouse")

# -----------------------------------------------------------------------------
# Part C: Combine and Finalize
# -----------------------------------------------------------------------------

# Combine both datasets
tf_targets <- bind_rows(human_data, mouse_data) %>%
  dplyr::rename(
    TF = tf,
    Target = target,
    Confidence = confidence,
    Mode_of_Regulation = mor
  ) %>%
  dplyr::select(Species, TF, Target, Confidence, Mode_of_Regulation) %>%
  as_tibble()

# -----------------------------------------------------------------------------
# Column Explanations:
# 1. Species: "Human" or "Mouse".
# 2. TF: Gene Symbol of the Transcription Factor.
# 3. Target: Gene Symbol of the Target Gene.
# 4. Confidence: Quality level (A=High, B=Moderate, C=Medium).
# 5. Mode_of_Regulation: 1 (Activation) or -1 (Inhibition).
# -----------------------------------------------------------------------------

# 4. Save the processed data to the package
# This creates 'data/tf_targets.rda'
usethis::use_data(tf_targets, overwrite = TRUE)

