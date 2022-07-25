# This script will be used to generate a table of the optimization results that is 
# sorted by log-likelihood and has a column for each parameter

# The script should be executed from the folder where all the optimization sub-folders
# are located.

## Load libraries: ##
library(tidyverse)

## Load functions: ##

# read an optimization table
read.opti.table <- function(file){
  table <- read_tsv(file, progress = F, show_col_types = FALSE)
  
  colnames(table)[7] <- "parameters"

  table$pop_pair <- unlist(strsplit(unlist(strsplit(file, "/"))[1], "_"))[5]
  
  table$opti_n <- unlist(strsplit(unlist(strsplit(file, "/"))[1], "_"))[6]
  
  return(table)
}

# format one row of the optimization table
format.opti.table <- function(df.row){
  raw <- df.row
  
  if (raw[1,1] == "model_1_a" | raw[1,1] == "model_1_b") {
    param_list <- c("Tsplit", "Tbot1", "iber_a", "iber_pr", "eura_a", "eura_pr", "m", "m_12", "m_21")

    raw$Tbot2 <- 0
    raw$iber_pr_a <- 1
    
  } else if (raw[1,1] == "model_2_a" | raw[1,1] == "model_2_b" | raw[1,1] == "model_2_c") {
    param_list <- c("Tsplit", "Tbot2", "Tbot1", "iber_a", "iber_pr_a", "iber_pr", "eura_a", "eura_pr", "m", "ma_12", "ma_21", "m_12", "m_21")
  }
  
  table <- raw %>%
    separate(., col=7, sep = ",", into = param_list)
  
  param_list_short <- c("Tsplit", "Tbot2", "Tbot1", "iber_a", "iber_pr_a", "iber_pr", "eura_a", "eura_pr")
  
  col_order <- c("Model", "pop_pair", "log-likelihood", "theta", param_list_short)
  
  table <- table[, col_order]
  
  return(table)
}

##

# list of optimization files
opti_files <- system("ls */*.optimized.txt", intern=T)

# empty df to be filled with optimization results
all_opti_likelihood <- data.frame()

# fill df with optimization results
for (opti_file in opti_files){
  all_opti_likelihood <- rbind(all_opti_likelihood, read.opti.table(opti_file))
}

# empty df to be filled with formatted optimization results
formatted_opti_table <- data.frame()

# fill df with formatted optimization results
for (n in 1:NROW(all_opti_likelihood)){
  row <- all_opti_likelihood[n,]
  formatted_opti_table <- rbind(formatted_opti_table, format.opti.table(row))
}

# order the formatted table based on log-likelihood
ordered_opti_table <- formatted_opti_table[order(-formatted_opti_table$`log-likelihood`),]

# save ordered and formatted table
write.table(ordered_opti_table, "ordered_opti_table.tsv", sep = '\t',
            row.names = F, quote = F)
## 
