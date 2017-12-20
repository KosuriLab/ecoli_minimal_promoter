# install.packages("plyr")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("cowplot")
# install.packages("lazyeval")
# install.packages("wesanderson")
# install.packages("tidyr")
# install.packages('stringr')

library("ggplot2")
library("dplyr")
library("wesanderson")
require(cowplot)
library("stringr")

options(stringsAsFactors = F)

split_exp <- read.table('../processed_data/rlp5Min_SplitVariants.txt', header = T)

#Figure 5A 
#Split data by -10 element and color by -35

#Organize data
pal <- wes_palette("Zissou", 8, type = "continuous")
legend_ord <- levels(with(split_exp, reorder(Minus35,-RNA_exp_average, median)))
split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))
split_exp$Minus10  <- with(split_exp, reorder(Minus10,RNA_exp_average, median))

ggplot(split_exp, aes(x = Minus10, y=RNA_exp_average)) +
  geom_jitter(aes(colour=Minus35), alpha = .5, size = 2) +
  geom_boxplot(alpha=0, lwd = .6) +
  scale_color_manual(values = pal, breaks = legend_ord) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10(limits = c(.05,39)) +
  labs(title = 'Distribution of Expression by -10', 
       y = 'Expression', x = '-10 Variant') +
  annotation_logticks(sides = 'l')
ggsave('../figs/Figure5A.pdf')

#Figure 5B

pal <- wes_palette("Zissou",10, type = "continuous")

split_exp$Minus10  <- with(split_exp, reorder(Minus10,RNA_exp_average, median))
split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))

split_exp  %>% 
    select(RNA_exp_average, Minus35, Minus10) %>% 
    filter(RNA_exp_average > 0) %>%
    group_by_(.dots=c("Minus35", "Minus10")) %>% 
    summarize(RNA_exp_average = median(RNA_exp_average)) %>% 
    ungroup() %>%
    ggplot(., aes(x=Minus10, y= Minus35)) + 
    geom_tile(aes(fill = log10(RNA_exp_average))) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = pal, name = 'Log10(Expression)', 
                         limits = c(-1.3,1.3)) +
    labs(title = 'Combinatorial Expression of Core Promoter Elements', 
         x = '-10 Variant', y = '-35 Variant')
ggsave('../figs/Figure5B.pdf')

#Figure5C

pal <- wes_palette("Zissou",8, type = "continuous")

split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(split_exp, reorder(Minus35,-RNA_exp_average, median)))

a <- split_exp %>% 
    filter(UP_element == 'gourse-326fold-up') %>% 
    select(RNA_exp_average, Minus35, Minus10, Spacer, Background) %>%
    filter(RNA_exp_average > 0) %>%
    ungroup()

b <-split_exp %>% 
    filter(UP_element == 'noUP') %>% 
    select(RNA_exp_average, Minus35, Minus10, Spacer, Background) %>%
    filter(RNA_exp_average > 0) %>%
    ungroup()

inner_join(a,b, by = c('Minus10', 'Minus35', 'Background', 'Spacer')) %>%
    transform(Fold_UP = RNA_exp_average.x/RNA_exp_average.y) %>% 
    ggplot(., aes(x=RNA_exp_average.y, y=Fold_UP)) + 
    geom_point(aes(color = Minus35)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
    ylim(0,25) +
    labs(title = 'Increase Due to UP-Element', 
         x = 'Base Promoter Strength', 
         y = 'Increase in Expression Due to UP Element')
ggsave("../figs/Figure5C.pdf")



#Figure 5D
#Look at Fold change due to UP-element for each combination of -10, -35

pal <- wes_palette("Zissou", 8, type = "continuous")
split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))
split_exp$Minus10  <- with(split_exp, reorder(Minus10,RNA_exp_average, median))

a<- split_exp %>% filter(UP_element == 'gourse-326fold-up') %>%
  select(RNA_exp_average, Minus35, Minus10) %>%
  filter(RNA_exp_average > 0) %>%
  group_by_(.dots=c("Minus35", "Minus10")) %>% 
  summarize(RNA_exp_average = median(RNA_exp_average)) %>%
  ungroup()

b <-split_exp %>% filter(UP_element == 'noUP') %>%
  select(RNA_exp_average, Minus35, Minus10) %>%
  filter(RNA_exp_average > 0) %>%
  group_by_(.dots=c("Minus35", "Minus10")) %>% 
  summarize(RNA_exp_average = median(RNA_exp_average)) %>%
  ungroup() 

inner_join(a,b, by = c('Minus10', 'Minus35')) %>%
  transform(Fold_UP = RNA_exp_average.x/RNA_exp_average.y) %>% 
  ggplot(., aes(x=Minus10, y= Minus35)) +
  geom_tile(aes(fill = log2(Fold_UP))) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = pal, name = 'Fold Increase') +
  labs(title = 'Fold Increase Due to UP-Element', x = '-10 Variant', y = '-35 Variant')
ggsave('../figs/Figure5D.pdf')





