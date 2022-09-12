# load libraries
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

opti_table_file = args[1]

# read optimization table
opti_table = read_tsv(opti_table_file)

# prepare dfs of all rounds
full_df = opti_table %>%
  separate(., col="Replicate", sep = "_", into = c(NA, "round", NA, "replicate"))

# colors for plotting
cols = brewer.pal(11, "Spectral")[c(1,4,9,11)]
full_df$color = ifelse(full_df$round == 1, cols[1], ifelse(full_df$round == 2, cols[2], ifelse(full_df$round == 3, cols[3], cols[4])))

# divide by model:
for (model in c('model_1_a', 'model_1_b', 'model_2_a', 'model_2_b', 'model_2_c')){
  print(model)

  # model only all rounds
  full_df_model = full_df %>% filter(Model == model)
  
  # model only and only 4th rounds
  fourth_round_df = full_df_model %>% filter(round == 4) 
  
  pops = unique(full_df_model$pop_pair)
    
  # dot plot with repetitions as facets
  # lolli_plot <- ggplot()+
  #   geom_point(data=full_df_model, aes(x=order, y=`log-likelihood`), size = 3,
  #              shape=21, fill = full_df_model$color) +
  #   facet_wrap(~ opti_n)
  # lolli_plot
# 
  
  # top 100 reconstructions ordered by likelihood
  full_df_model$ll_order = c(1:NROW(full_df_model))
  full_df_model$opti_n = as.factor(full_df_model$opti_n)
  
  barlolli_plot <- ggplot(head(full_df_model, 10), aes(x=-ll_order, y=`log-likelihood`, fill = opti_n)) + geom_bar(stat ="identity") +
    scale_fill_brewer(palette="Set1")
  
  ggsave(filename = paste0("plots/", pops, "-", model, "_bar_top10_lolli_plot.pdf"), plot = barlolli_plot,
         width = 8, height = 8)
  
  # line plot with repetitions as facets
  linylolli_plot <- ggplot()+
    geom_line(data=full_df_model, aes(x=order, y=`log-likelihood`), size = 1) +
    facet_wrap(~ opti_n)
  
  ggsave(filename = paste0("plots/", pops, "-", model, "_lolli_plot.pdf"), plot = linylolli_plot,
         width = 8, height = 8)
  
  lolli_plot_fourth <- ggplot()+
    geom_line(data=fourth_round_df, aes(x=order, y=`log-likelihood`), size = 1) +
    facet_wrap(~ opti_n)
  
  ggsave(filename = paste0("plots/", pops, "-", model, "_lolli_fourth_panels_plot.pdf"), plot = lolli_plot_fourth,
         width = 8, height = 8)
  
  # dot plot only fourth round all together
  lolli_plot_fourth_together <- ggplot()+
    geom_point(data=fourth_round_df, aes(x=opti_n, y=`log-likelihood`), size = 3,
               shape=21, fill = fourth_round_df$opti_n)
  
  ggsave(filename = paste0("plots/", pops, "-", model, "_lolli_fourth_together_plot.pdf"), plot = lolli_plot_fourth_together,
         width = 8, height = 8)
}
