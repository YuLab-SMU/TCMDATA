#' @description add microbes-metabolites-targets relationship data from gutMGene
#' @references Qi, C., He, G., Qian, K., Guan, S., Li, Z., Liang, S., Liu, J., Ke, X., Zhang, S., Lu, M., Cheng, L., & Zhang, X. (2025). 
#' gutMGene v2.0: an updated comprehensive database for target genes of gut microbes and microbial metabolites. Nucleic acids research, 53(D1), D783–D788. 
#' https://doi.org/10.1093/nar/gkae1002


rm(list = ls())
library(dplyr)
## load raw data
df_bac_met <- read.csv("/Users/hinna/Desktop/gutMGene/Gut Microbe-Microbial metabolite.csv")
df_bac_gene <- read.csv("/Users/hinna/Desktop/gutMGene/Gut Microbe-Host Gene.csv")
df_met_gene <- read.csv("/Users/hinna/Desktop/gutMGene/Microbial metabolite-Host Gene.csv")

## process function
process_species_data <- function(species_name) {
  
  # microbe-metablites
  clean_bac_met <- df_bac_met %>%
    filter(human.mouse == species_name) %>% 
    select(
      Bacteria = Gut.Microbiota,
      Bacteria_ID = Gut.Microbiota.NCBI.ID,
      Metabolite,
      Metabolite_ID = Metabolite.PubChem.CID,
      PMID_Bac_Met = PMID
    ) %>%
    distinct() 
  
  # metabolite-targets
  clean_met_gene <- df_met_gene %>%
    filter(human.mouse == species_name) %>% 
    select(
      Metabolite, 
      Target = Gene,
      Target_ID = Gene.ID,
      Interaction = Alteration, 
      PMID_Met_Target = PMID
    ) %>%
    distinct()
  
  merged_axis <- full_join(clean_bac_met, clean_met_gene, by = "Metabolite") %>%
    select(
      Bacteria, 
      Bacteria_ID,
      Metabolite, 
      Metabolite_ID, 
      Target, 
      Target_ID, 
      Interaction,
      PMID_Bac_Met,
      PMID_Met_Target)
  
  return(merged_axis)
}

## process raw data
gut_axis_human <- process_species_data("human")
gut_axis_mouse <- process_species_data("mouse")
gutMGene <- list(gut_axis_human = gut_axis_human,
                 gut_axis_mouse = gut_axis_mouse)

usethis::use_data(gutMGene, overwrite = TRUE)

