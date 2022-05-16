# load libraries
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(gridExtra)

# read input args
args = commandArgs(trailingOnly=TRUE)
# define models
model <- args[1]
# define population
pops <- args[2]

# create df for all model's replicates
full_df <- data.frame()

for(r in 1:20){
  print(r)
  filename <- paste0("~/Documents/Introgression/dadi/tables/", pops, ".",
                     model, ".optimized.r_", r, ".txt")
  if(file.exists(filename) == "FALSE") {
    next
  }
  opti_table <- read_tsv(filename) %>%
    separate(., col="Replicate", sep = "_", into = c(NA, "round", NA, "replicate"))
  
  cols <- brewer.pal(11, "Spectral")[c(1,4,9,11)]
  
  opti_table$color <- ifelse(opti_table$round == 1, cols[1], ifelse(opti_table$round == 2, cols[2], ifelse(opti_table$round == 3, cols[3], cols[4])))
  opti_table$repetition <- r
  opti_table$replicate <- c(1:NROW(opti_table))
  
  full_df <- rbind(full_df, opti_table)
  
}

# only fourth round df
fourth_round_df <- full_df %>% filter(round == 4) 

# dot plot with repetitions as facets
lolli_plot <- ggplot()+
  geom_point(data=full_df, aes(x=replicate, y=`log-likelihood`), size = 3,
             shape=21, fill = full_df$color) +
  facet_wrap(~ repetition)

# line plot only fourth round facets
lolli_plot_fourth <- ggplot()+
  geom_line(data=fourth_round_df, aes(x=replicate, y=`log-likelihood`), size = 1) +
  facet_wrap(~ repetition)

# dot plot only fourth round all together
lolli_plot_fourth_together <- ggplot()+
  geom_point(data=fourth_round_df, aes(x=replicate, y=`log-likelihood`), size = 3,
             shape=21, fill = fourth_round_df$repetition)

pdf(file = paste0("plots/", pops, "-", model, "_lolli_plot.pdf"),
    width = 8,
    height = 8)
lolli_plot
dev.off()

pdf(file = paste0("plots/", pops, "-", model, "_lolli_fourth_panels_plot.pdf"),
    width = 8,
    height = 8)
lolli_plot_fourth
dev.off()

pdf(file = paste0("plots/", pops, "-", model, "_lolli_fourth_together_plot.pdf"),
    width = 8,
    height = 8)
lolli_plot_fourth_together
dev.off()
