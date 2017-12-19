# install.packages("plyr")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("cowplot")
# install.packages("lazyeval")
# install.packages("wesanderson")
# install.packages("tidyr")
#install.packages('stringr')

library("ggplot2")
library("dplyr")
library("wesanderson")
require(cowplot)
library("stringr")

options(stringsAsFactors = F)

split_exp <- read.table('./rlp5Min_SplitVariants.txt', header = T)

#Figure 6A
# Split data by background, look at distribution of expression of each background while highlighting the -10

split_exp$Background  <- with((split_exp %>% mutate(GC_Content = (((str_count(variant, pattern = "G")) + (str_count(variant, pattern = "C")))/150))), reorder(Background,GC_Content, median))

ggplot(split_exp, aes(x= Background, y=RNA_exp_average)) + stat_boxplot(geom = "errorbar", width = 0.25) + geom_boxplot(lwd = .6, fill = "#3B9AB2") +  
  geom_jitter(data=filter(split_exp, Minus10 == 'consensus10' & Minus35 == 'consensus35'), aes(x = Background, y=RNA_exp_average), color = "#F21A00", alpha=.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + annotation_logticks(sides = 'y') + scale_y_log10(limits = c(.05,39)) +
  labs(title = 'Background Expression profiles with Consensus -10,-35', y = 'Expression', x = 'Background Variant') + annotation_logticks(sides = 'l') 
ggsave('Figure6A.pdf')

#Figure 6B
#Here we take all promoters with expression > .5 (This is expression above the negative controls), 
#calculate the GC content of the 8bp of the spacer proximal to the -10,  
#and calculate correlation between median RNA expression of each spacer and GC content 

spacerseqs <- read.table("./SpacerSequences.txt", header = T, sep = "\t")

spacerseqs <- spacerseqs %>% mutate(endSeqs = (((str_count(substr(Sequence,10,17), pattern = "G")) + (str_count(substr(Sequence,10,17), pattern = "C")))/(17-9)), Length = (17-9)) 

temp <- split_exp %>% select(RNA_exp_average, Spacer) %>% 
  group_by(Spacer) %>%  
  filter(RNA_exp_average > .5) %>%
  summarize(RNA_med_1 = median(RNA_exp_average)) %>%
  left_join(spacerseqs, by = 'Spacer') %>% ungroup()

#calculate correlation
corr <- cor(temp$RNA_med_1, temp$endSeqs)
cor.test(temp$RNA_med_1, temp$endSeqs) #pvalue is .0358

split_exp %>% select(RNA_exp_average, Spacer) %>% group_by(Spacer) %>%  
  filter(RNA_exp_average > .5) %>%
  summarize(RNA_med_1 = median(RNA_exp_average)) %>%
  left_join(spacerseqs, by = 'Spacer') %>% 
  ungroup() %>%
  ggplot(aes(x=endSeqs, y=RNA_med_1)) +
  geom_point(aes(color = "#3B9AB2"),  size = 3) +
  geom_smooth(method=lm) + 
  annotate("text", x =.7, y = 4, label = paste('r==', signif(corr, 3)), parse = T) +
  geom_text(aes(label=Spacer),hjust=.05, vjust=-1) + 
  scale_color_manual(values = pal, name = 'Spacer') + 
  xlim(0,.85) +
  labs(x='%GC Content', y='Average Spacer Expression', title='GC Content of Spacer Correlates with Promoter Expression') 
ggsave('Figure6B.pdf')
