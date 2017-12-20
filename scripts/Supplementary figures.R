#This script contains the code used to generate Supplementary Figures

library("ggplot2")
library("dplyr")
library("wesanderson")
library(stringr)
library(nnet)
library(devtools)
require(cowplot)

options(stringsAsFactors = F)

split_exp <- read.table('../processed_data/rlp5Min_SplitVariants.txt', header = T)
min_MOPS <- read.table("../processed_data/revLP5_Min_MOPS_glu_expression.txt", header = T)

pal <- wes_palette("Zissou",8, type = "continuous")

#Supplementary figure 2

ggplot(min_MOPS, aes(DNA_sum)) + 
    geom_density(alpha=.2, fill = "#3B9AB2", color = "#3B9AB2") +
    scale_x_log10() + annotation_logticks(sides = 'b') +
    annotate("text", x=10, y=1.35, 
             label = paste('Median Normalized RPM per variant =', 
                           signif(median(min_MOPS$DNA_sum), 2))) +
    labs(x='DNA RPM per Variant', y='Density', 
         title = 'Distribution of the Number of Reads Per Variant')
ggsave("../figs/SupplementaryFigure2.png")


#Supplementary Figure 3
#This figure is a combination of four graphs 

split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))
split_exp$Minus10 <- with(split_exp, reorder(Minus10,-RNA_exp_average, median))
split_exp$Background <- with(split_exp, reorder(Background,RNA_exp_average, median))
split_exp$Spacer <- with(split_exp, reorder(Spacer,RNA_exp_average, median))

#NO up
ggplot(filter(split_exp, UP_element == 'noUP'), aes(x = Spacer, y = Background)) + 
  geom_tile(aes(fill = log10(RNA_exp_average))) + 
  scale_fill_gradientn(colors = pal, limits = c(-1.4,1.6)) + 
  facet_grid(Minus10 ~ Minus35) +
  theme(panel.background = element_rect(fill = "gray80"), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        strip.text = element_blank(),
        panel.spacing = unit(0.10, "lines"))
ggsave('../figs/SF3_alldata_noup.pdf')

#136x UP
ggplot(filter(split_exp, UP_element == 'gourse-136fold-up'), aes(x = Spacer, y = Background)) + 
  geom_tile(aes(fill = log10(RNA_exp_average))) + 
  scale_fill_gradientn(colors = pal, limits = c(-1.4,1.6)) + 
  facet_grid(Minus10 ~ Minus35) +
  theme(panel.background = element_rect(fill = "gray80"), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        strip.text = element_blank(),
        panel.spacing = unit(0.10, "lines"))
ggsave('../figs/SF3_alldata_136.pdf')

#326x UP
ggplot(filter(split_exp, UP_element == 'gourse-326fold-up'), aes(x = Spacer, y = Background)) + 
  geom_tile(aes(fill = log10(RNA_exp_average))) + 
  scale_fill_gradientn(colors = pal, limits = c(-1.4,1.6)) + 
  facet_grid(Minus10 ~ Minus35) +
  theme(panel.background = element_rect(fill = "gray80"), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        strip.text = element_blank(),
        panel.spacing = unit(0.10, "lines"))
ggsave('../figs/SF3_alldata_326.pdf')


#Single example
ggplot(filter(split_exp, Minus10 == 'consensus10', Minus35 == 'consensus35', UP_element == 'noUP'), aes(x = Spacer, y = Background)) + 
  geom_tile(aes(fill = log10(RNA_exp_average))) + 
  scale_fill_gradientn(colors = pal, limits = c(-1.4,1.6)) + 
  theme(panel.background = element_rect(fill = "gray80"), 
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 6), 
        axis.text.x = element_text(size = 6, angle = -90), 
        strip.text = element_text(size = 6),
        panel.spacing = unit(0.40, "lines"))
ggsave('../figs/SF3_noUPdata.pdf')


#Supplementary Figure 4

split_exp <- mutate(split_exp, RNA_exp_log = log10(RNA_exp_average))

#Sample training and test data

#split_exp <- mutate(split_exp, RNA_exp_average = log10(RNA_exp_average)) #Log transformed data for lm

smp_size <- floor(0.5 * nrow(split_exp))

set.seed(123)
train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)

train <- split_exp[train_ind, ]
test <- split_exp[-train_ind, ]

#supp4A, NO log transformation, no Interaction term
linM = train %>%
  select(RNA_exp_average, UP_element, Minus35, Spacer, Minus10, Background) %>%
  lm(formula=RNA_exp_average ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=., model = T) 

test$predicted_exp_lm <- predict(linM, test) #predict expression of test data

corr <- summary(lm(RNA_exp_average ~ predicted_exp_lm, test))$r.square

test$Minus35  <- with(test, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_average, median)))


ggplot(test, aes(predicted_exp_lm, RNA_exp_average, color = Minus35)) +
  geom_point(alpha = .4, aes(color = Minus35)) +
  annotate("text", x = .1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
  scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
  scale_x_log10(limits = c(0.05,30)) + scale_y_log10(limits = c(0.05,30)) + annotation_logticks() + 
  xlab('Predicted Expression') + ylab('Expression') +
  ggtitle('no log linear, no interaction term model') 
ggsave('../figs/SF4_nolog-linear, no interaction.pdf')


#supp4B,Linear model with interaction term
linM = train %>%
  select(RNA_exp_average, UP_element, Minus35, Spacer, Minus10, Background) %>%
  lm(formula=RNA_exp_average ~ UP_element + Minus35 + Spacer + Minus10 + Background + Minus10:Minus35, 
     data=., model = T) 

test$predicted_exp_lm <- predict(linM, test) #predict expression of test data

corr <- summary(lm(RNA_exp_average ~ predicted_exp_lm, test))$r.square

test$Minus35  <- with(test, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_average, median)))

ggplot(test, aes(predicted_exp_lm, RNA_exp_average, color = Minus35)) +
  geom_point(alpha = .4, aes(color = Minus35)) +
  annotate("text", x = .1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
  scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
  scale_x_log10(limits = c(0.05,30)) + scale_y_log10(limits = c(0.05,30)) + annotation_logticks() + 
  xlab('Predicted Expression') + ylab('Expression') +
  ggtitle('no log linear model') 
ggsave('../figs/SF4_nolog-linear.pdf')

#Supp4C, log-linear modelWithout Interaction term
linM = train %>%
  select(RNA_exp_log, UP_element, Minus35, Spacer, Minus10, Background) %>%
  lm(formula=RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=., model = T) 

test$predicted_exp_lm <- predict(linM, test) #predict expression of test data

corr <- summary(lm(RNA_exp_log ~ predicted_exp_lm, test))$r.square

test$Minus35  <- with(test, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_average, median)))

ggplot(test, aes(10^predicted_exp_lm, 10^RNA_exp_log, color = Minus35)) +
    geom_point(alpha = .4, aes(color = Minus35)) +
    annotate("text", x = .1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
    scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
    scale_x_log10(limits = c(0.05,30)) + scale_y_log10(limits = c(0.05,30)) + 
    annotation_logticks() + 
    xlab('Predicted Expression') + ylab('Expression') +
    ggtitle('log-linear no interaction term') 
ggsave('../figs/SF4_log-linear_noint.pdf')


#Supp4D, neural network with only 5% of training data

smp_size <- floor(0.05 * nrow(split_exp))

set.seed(124)
train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)

train <- split_exp[train_ind, ]
test <- split_exp[-train_ind, ]

set.seed(124)
fit <- nnet(RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background,
            data=train, size=10, maxit=300, linout=T, decay=0.01)

# make predictions
test$predicted_exp_NN <- predict(fit, test, type="raw")

corr <- summary(lm(RNA_exp_log ~ predicted_exp_NN, test))$r.square

test$Minus35  <- with(test, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_average, median)))

ggplot(test, aes(10^predicted_exp_NN, 10^RNA_exp_log, color = Minus35)) + 
  geom_point(alpha = .4) +
  annotate("text", x =0.1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
  scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) + 
  scale_x_log10(limits = c(0.05,30)) + scale_y_log10(limits = c(0.05,30)) + annotation_logticks() + 
  xlab('Predicted Expression') +
  ylab('Expression') +
  ggtitle('Neural Network on .05')
ggsave('../figs/SF4_5percent_NN.pdf')

#Supplementary Figure 5, Neural Network schematic
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
plot.nnet(fit)


#Supplementary Figure 6
#THis figure is in two parts
split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(split_exp, reorder(Minus35,-RNA_exp_average, median)))

split_exp$Minus10 <- with(split_exp, reorder(Minus10,RNA_exp_average, median))
legend_ord <- levels(with(split_exp, reorder(Minus10,-RNA_exp_average, median)))

#326xUP Log2 (TOP)
split_exp %>% 
    filter(UP_element == 'gourse-326fold-up') %>% 
    select(RNA_exp_average, Minus35, Minus10) %>% 
    filter(RNA_exp_average > 0) %>%
    group_by_(.dots=c("Minus35", "Minus10")) %>% 
    summarize(RNA_exp_average = median(RNA_exp_average)) %>% ungroup() %>%
    ggplot(., aes(x=Minus10, y= Minus35)) +
    geom_tile(aes(fill = log2(RNA_exp_average))) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = pal, name = 'Log2(Expression)', limits = c(-4,4)) +
    labs(title = 'Combinatorial Expression of Core Promoter Elements (326xUP Element)', 
         x = '-10 Variant', y = '-35 Variant')
ggsave("../figs/SF6_Combinatorial_326xUP_log2.pdf")

#NO UP (Bottom)
split_exp %>% 
    filter(UP_element == 'noUP') %>% 
    select(RNA_exp_average, Minus35, Minus10) %>% 
    filter(RNA_exp_average > 0) %>%
    group_by_(.dots=c("Minus35", "Minus10")) %>% 
    summarize(RNA_exp_average = median(RNA_exp_average)) %>% ungroup() %>%
    ggplot(., aes(x=Minus10, y= Minus35)) +
    geom_tile(aes(fill = log2(RNA_exp_average))) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = pal, name = 'Log2(Expression)', limits = c(-4,4)) +
    labs(title = 'Combinatorial Expression of Core Promoter Elements (326xUP Element)', 
         x = '-10 Variant', y = '-35 Variant')
ggsave("../figs/SF6_Combinatorial_noUP_log2.pdf")

#Supplementary Figure 7 

split_exp$Minus10 <- with(split_exp, reorder(Minus10,RNA_exp_average, median))
legend_ord <- levels(with(split_exp, reorder(Minus10,-RNA_exp_average, median)))

a <- split_exp %>% 
    filter(UP_element == 'gourse-326fold-up') %>% 
    select(RNA_exp_average, Minus35, Minus10, Spacer, Background) %>%
    filter(RNA_exp_average > -0) %>%
    ungroup()

b <-split_exp %>% 
    filter(UP_element == 'noUP') %>% 
    select(RNA_exp_average, Minus35, Minus10, Spacer, Background) %>%
    filter(RNA_exp_average > -0) %>%
    ungroup()

inner_join(a,b, by = c('Minus10', 'Minus35', 'Background', 'Spacer')) %>%
  transform(Fold_UP = RNA_exp_average.x/RNA_exp_average.y) %>% 
  ggplot(., aes(x=RNA_exp_average.y, y=Fold_UP)) + geom_point(aes(color = Minus10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + annotation_logticks(sides = 'l') +
  scale_color_manual(values = pal, name = '-10 Variant', breaks = legend_ord) +
  ylim(0,25) + scale_y_log10() + geom_hline(yintercept=1, color = 'black') +
  labs(title = 'Increase Due to UP-Element', x = 'Base Promoter Strength', 
       y = 'Increase in Expression Due to UP Element')
ggsave("../figs/SF7_UPIncrease_log.pdf")


#Supplementary Figure 8


split_exp$Background <- with((split_exp %>% mutate(GC_Content = (((str_count(variant, pattern = "G")) + (str_count(variant, pattern = "C")))/150))), reorder(Background,GC_Content, median))

split_exp %>% select(RNA_exp_average, Minus35, Minus10, Background) %>% filter(RNA_exp_average > 0) %>%
  group_by_(.dots=c("Minus35", "Minus10", "Background")) %>% 
  summarize(RNA_exp_average = median(RNA_exp_average)) %>% ungroup() %>%
  ggplot(., aes(x=Minus10, y= Minus35)) + 
  geom_tile(aes(fill = log10(RNA_exp_average))) +
  #ggplot(filter(., Background != 'bg4323949:4324099', Background != 'bg977040:977190'), aes(x = Spacer, y = Background)) + 
  #geom_tile(aes(fill = log10(RNA_exp_average))) + 
  scale_fill_gradientn(colors = pal, limits = c(-1.3,1.3)) + 
  facet_wrap(~Background, ncol = 4) +
  theme(panel.background = element_rect(fill = "gray80"), 
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 6), 
        axis.text.x = element_text(size = 6, angle = -90), 
        strip.text = element_text(size = 10),
        panel.spacing = unit(0.20, "lines"))
ggsave("../figs/PreferredCorePromoterperBG.pdf")

#Supplementary Figure 9

split_exp <- mutate(split_exp, RNA_exp_log = log10(RNA_exp_average)) #Log transformed data for lm

smp_size <- floor(0.5 * nrow(split_exp))

set.seed(123)
train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)

train <- split_exp[train_ind, ]
test <- split_exp[-train_ind, ]

set.seed(123)
fit <- nnet(RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10, 
            data=train, size=10, maxit=300, linout=T, decay=0.01)

# make predictions
test$predicted_exp_NN <- predict(fit, test, type="raw")

corr <- summary(lm(RNA_exp_log ~ predicted_exp_NN, test))$r.square
summary(lm(RNA_exp_log ~ predicted_exp_NN, test)) #p < 2.2 x 1016

test$Background  <- with(test, reorder(Background,RNA_exp_average, median))
legend_ord <- levels(with(test, reorder(Background,-RNA_exp_average, median)))

ggplot(test, aes(10^predicted_exp_NN, 10^RNA_exp_log, color = Background)) + 
  geom_point(alpha = .4) +
  annotate("text", x =0.1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T, size = 6) +
  scale_color_manual(values = pal, name = 'Background', breaks = legend_ord) + 
  scale_x_log10(limits = c(0.05,30)) + scale_y_log10(limits = c(0.05,30)) + annotation_logticks() + 
  xlab('Predicted Expression') +
  ylab('Expression') +
  ggtitle('Using Neural Network to Predict Variant Expression')
ggsave('../figs/SF9_NN_noBG.png')


#Supplementary Figure 10

split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(split_exp, reorder(Minus35,-RNA_exp_average, median)))

split_exp$Minus10 <- with(split_exp, reorder(Minus10,RNA_exp_average, median))
legend_ord <- levels(with(split_exp, reorder(Minus10,-RNA_exp_average, median)))

#ECK726, TOP
split_exp %>% filter(Spacer == 'ECK125137726') %>%
  select(RNA_exp_average, Minus10, Minus35) %>%
  filter(RNA_exp_average > 0) %>%
  group_by_(.dots=c("Minus10", "Minus35")) %>% 
  summarize(RNA_exp_average = median(RNA_exp_average)) %>%
  ungroup() %>%
  ggplot(., aes(x=Minus10, y=Minus35)) +
  geom_tile(aes(fill = log10(RNA_exp_average))) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = pal, name = 'Log10(Expression)', limits = c(-1.3,1.3)) +
  labs(title = 'ECK726', x = '-10 Variant', y = '-35 Variant')
ggsave('../figs/SF10_ECK726.pdf')

#ECK405, Bottom
split_exp %>% filter(Spacer == 'ECK125137405') %>%
  select(RNA_exp_average, Minus10, Minus35) %>%
  filter(RNA_exp_average > 0) %>%
  group_by_(.dots=c("Minus10", "Minus35")) %>% 
  summarize(RNA_exp_average = median(RNA_exp_average)) %>%
  ungroup() %>% ggplot(., aes(x=Minus10, y=Minus35)) +
  geom_tile(aes(fill = log10(RNA_exp_average))) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = pal, name = 'Log10(Expression)', limits = c(-1.3,1.3)) +
  labs(title = 'ECK405', x = '-10 Variant', y = '-35 Variant')
ggsave('../figs/SF10_ECK405.pdf')


# Supplementary Figure X, Kinney energy matrix
min_RNAP_Scores <- read.table('../processed_data/min_rnap_scores.txt', 
                              sep = '\t', header = T)
#Get correlation
corr <- summary(lm(rnap_site_no_spacer_score ~ log10(RNA_exp_average), min_RNAP_Scores))$r.squared

#Order data so -35 color is consistent
legend_ord <- levels(with(min_RNAP_Scores, reorder(Minus35,-RNA_exp_average, median)))

#Plot, colored by -35
ggplot(min_RNAP_Scores, aes(rnap_site_no_spacer_score, RNA_exp_average, color = Minus35)) +
    geom_point(alpha = .4, aes(color = Minus35)) +
    annotate("text", x = 110, y = 1.5, label = paste('R^2==', signif(corr, 3)), parse = T) +
    scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
    annotation_logticks(sides = 'l') + scale_y_log10(limits = c(.03,30)) +
    xlab('RNAP Site Score') + ylab('Expression') +# coord_flip(xlim=c(130,30)) + 
    scale_x_reverse(lim=c(130,0)) +
    ggtitle('Comparing RNAP Binding Energy to Expression') 
ggsave('../figs/SFX_kinney_matrix.pdf')
