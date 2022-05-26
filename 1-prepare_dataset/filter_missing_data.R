# load libraries
library(tidyverse)
library(RColorBrewer)

# input file from pop_vcf_missing_gts.sh
filename <- "1-prepare_dataset/tables/lp_ll_introgression_perpop_missing_gts.csv"

# load table
miss_table <- read_csv(filename, progress = T, show_col_types = FALSE)

# extract each population's counts of missing values and calculate missing proportion
# by dividing by number of missing samples (Var1) by total number of samples in population:
freq_miss_table <- data.frame()

for (pop in c("lpa", "wel", "eel", "sel")) {
  print(pop)  
  pop_miss_table <- data.frame(table(miss_table[,pop]))
  pop_miss_table$Var1 <- as.numeric(levels(pop_miss_table$Var1))/max(as.numeric(levels(pop_miss_table$Var1)))
  pop_miss_table$Freq_prop <- (pop_miss_table$Freq)/NROW(miss_table)
  pop_miss_table$Freq_cumsum <- cumsum(pop_miss_table$Freq_prop)
  pop_miss_table$population <- pop
  freq_miss_table <- rbind(freq_miss_table, pop_miss_table)
}

# extract each population into its own table for easier plot
lpa <- freq_miss_table %>% filter(population == "lpa")
wel <- freq_miss_table %>% filter(population == "wel")
eel <- freq_miss_table %>% filter(population == "eel")
sel <- freq_miss_table %>% filter(population == "sel")

# cumulative proportion of SNPs included vs decreasing missingness filter strictness
cummiss <- ggplot() +
  geom_line(data=lpa, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="blue", size=1.5)+
  geom_point(data=lpa, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="blue", size=4, shape=15)+
  geom_line(data=wel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="red", size=1.5)+
  geom_point(data=wel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="red", size=4, shape=15)+
  geom_line(data=eel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="black", size=1.5)+
  geom_point(data=eel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="black", size=4, shape=15)+
  geom_line(data=sel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="green", size=1.5)+
  geom_point(data=sel, aes(x=Var1, y=Freq_cumsum), alpha=0.5, color="green", size=4, shape=15)+
  scale_x_continuous(limits = c(0, 0.4))+
  xlab("Proportion of missing data")+
  ylab("Proportion of SNPs included")+
  theme_minimal()

write.table(x = freq_miss_table,
            file = paste0("1-prepare_dataset/tables/freq_miss_table.tsv"),
            quote=FALSE,  col.names = T, row.names = FALSE, sep= "\t")

pdf(file = paste0("1-prepare_dataset/plots/cumulative_miss.pdf"), width = 8, height = 4)
cummiss
dev.off()
