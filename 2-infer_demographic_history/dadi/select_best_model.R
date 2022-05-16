library(tidyverse)
library(RColorBrewer)

models <- c("model_1_a", "model_1_b", "model_2_a", "model_2_b", "model_2_c", "model_2_d")

pops <- "sm_ki"

full_df <- data.frame()

for (model in models){
  
  for(r in 1:20){

    filename <- paste0("~/Documents/Introgression/dadi/tables/", pops, ".",
                     model, ".optimized.r_", r, ".txt")
    if(file.exists(filename) == "FALSE") {
      next
    }
    opti_table <- read_tsv(filename, progress = F, show_col_types = FALSE) %>%
      separate(., col="Replicate", sep = "_", into = c(NA, "round", NA, "replicate"))
    colnames(opti_table)[8] <- "optimized_params"
    cols <- brewer.pal(11, "Spectral")[c(1,4,9,11)]
  
    opti_table$color <- ifelse(opti_table$round == 1, cols[1], ifelse(opti_table$round == 2, cols[2], ifelse(opti_table$round == 3, cols[3], cols[4])))
    opti_table$repetition <- r
    opti_table$replicate <- c(1:NROW(opti_table))
  
    full_df <- rbind(full_df, opti_table)
  
  }
}

full_df_ordered <- full_df[order(-full_df$`log-likelihood`),]

best_scoring_ll_row <- full_df_ordered[c(1),]

# TRANSLATE PARAMETERS:
