# This script contains the code used to generate Supplementary Figures
library("ggplot2")
library("dplyr")
library("wesanderson")
require(cowplot)

options(stringsAsFactors = F)

split_exp <- read.table('./rlp5Min_SplitVariants.txt', header = T)
min_MOPS <- read.table("./revLP5_Min_MOPS_glu_expression.txt", header = T)

# Supplementary figure 2
ggplot(min_MOPS, aes(DNA_sum)) + 
    geom_density(alpha=.2, fill = "#3B9AB2", color = "#3B9AB2") +
    scale_x_log10() + annotation_logticks(sides = 'b') +
    annotate("text", x=4, y=1.35, label = paste('Median Normalized RPM per variant =', signif(median(min_MOPS$DNA_sum), 2))) +
    labs(x='DNA RPM per Variant', y='Density', title = 'Distribution of the Number of Reads Per Variant')
ggsave("SupplementaryFigure2.png")

# Supplementary Figure 3
# This figure is a combination of four graphs 
split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_12, median))
split_exp$Minus10 <- with(split_exp, reorder(Minus10,-RNA_exp_12, median))
split_exp$Background <- with(split_exp, reorder(Background,RNA_exp_12, median))
split_exp$Spacer <- with(split_exp, reorder(Spacer,RNA_exp_12, median))

# NO up
ggplot(filter(split_exp, UP_element == 'noUP'), aes(x = Spacer, y = Background)) + 
    geom_tile(aes(fill = log10(RNA_exp_12))) + 
    scale_fill_gradientn(colors = pal, limits = c(-1.4,1.3)) + 
    facet_grid(Minus10 ~ Minus35) +
    theme(panel.background = element_rect(fill = "gray80"), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(), 
          axis.text.x = element_blank(), 
          strip.text = element_blank(),
          panel.spacing = unit(0.10, "lines"))
ggsave('alldata_noup.pdf')

# 136x UP
ggplot(filter(split_exp, UP_element == 'gourse-136fold-up'), aes(x = Spacer, y = Background)) + 
    geom_tile(aes(fill = log10(RNA_exp_12))) + 
    scale_fill_gradientn(colors = pal, limits = c(-1.4,1.3)) + 
    facet_grid(Minus10 ~ Minus35) +
    theme(panel.background = element_rect(fill = "gray80"), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(), 
          axis.text.x = element_blank(), 
          strip.text = element_blank(),
          panel.spacing = unit(0.10, "lines"))
ggsave('alldata_136.pdf')

# 326x UP
ggplot(filter(split_exp, UP_element == 'gourse-326fold-up'), aes(x = Spacer, y = Background)) + 
    geom_tile(aes(fill = log10(RNA_exp_12))) + 
    scale_fill_gradientn(colors = pal, limits = c(-1.4,1.3)) + 
    facet_grid(Minus10 ~ Minus35) +
    theme(panel.background = element_rect(fill = "gray80"), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(), 
          axis.text.x = element_blank(), 
          strip.text = element_blank(),
          panel.spacing = unit(0.10, "lines"))
ggsave('alldata_326.pdf')

# Single example
ggplot(filter(split_exp, Minus10 == 'consensus10', Minus35 == 'consensus35', UP_element == 'noUP'), aes(x = Spacer, y = Background)) + 
    geom_tile(aes(fill = log10(RNA_exp_12))) + 
    scale_fill_gradientn(colors = pal, limits = c(-1.4,1.3)) + 
    theme(panel.background = element_rect(fill = "gray80"), 
          axis.ticks = element_blank(), 
          axis.text.y = element_text(size = 6), 
          axis.text.x = element_text(size = 6, angle = -90), 
          strip.text = element_text(size = 6),
          panel.spacing = unit(0.40, "lines"))
ggsave('noUPdata.pdf')

# Supplementary Figure 4
split_exp <- mutate(split_exp, RNA_exp_log = log10(RNA_exp_12))

# Sample training and test data
# split_exp <- mutate(split_exp, RNA_exp_12 = log10(RNA_exp_1)) #Log transformed data for lm
smp_size <- floor(0.5 * nrow(split_exp))
set.seed(123)
train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)

train <- split_exp[train_ind, ]
test <- split_exp[-train_ind, ]

# supp4A, NO log transformation, no Interaction term
linM = train %>%
    select(RNA_exp_12, UP_element, Minus35, Spacer, Minus10, Background) %>%
    lm(formula=RNA_exp_12 ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=., model = T) 

test$predicted_exp_lm <- predict(linM, test) #predict expression of test data

corr <- summary(lm(RNA_exp_12 ~ predicted_exp_lm, test))$r.square

pal <- wes_palette("Zissou", 8, type = "continuous")
test$Minus35  <- with(test, reorder(Minus35,RNA_exp_1, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_1, median)))


ggplot(test, aes(predicted_exp_lm, RNA_exp_12, color = Minus35)) +
    geom_point(alpha = .4, aes(color = Minus35)) +
    annotate("text", x = .05, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
    scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
    scale_x_log10(limits = c(0.01,30)) + scale_y_log10(limits = c(0.01,30)) + 
    annotation_logticks() + 
    xlab('Predicted Expression') + ylab('Expression') +
    ggtitle('no log linear, no interaction term model') 
ggsave('nolog-linear, no interaction.pdf')

# supp4B,Linear model with interaction term
linM = train %>%
    select(RNA_exp_12, UP_element, Minus35, Spacer, Minus10, Background) %>%
    lm(formula=RNA_exp_12 ~ UP_element + Minus35 + Spacer + Minus10 + Background + Minus10:Minus35, data=., model = T) 

test$predicted_exp_lm <- predict(linM, test) #predict expression of test data

corr <- summary(lm(RNA_exp_12 ~ predicted_exp_lm, test))$r.square

pal <- wes_palette("Zissou", 8, type = "continuous")
test$Minus35  <- with(test, reorder(Minus35,RNA_exp_1, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_1, median)))

ggplot(test, aes(predicted_exp_lm, RNA_exp_12, color = Minus35)) +
    geom_point(alpha = .4, aes(color = Minus35)) +
    annotate("text", x = .05, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
    scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
    scale_x_log10(limits = c(0.01,30)) + scale_y_log10(limits = c(0.01,30)) + 
    annotation_logticks() + 
    xlab('Predicted Expression') + ylab('Expression') +
    ggtitle('no log linear model') 
ggsave('nolog-linear.pdf')

#Supp4C, log-linear modelWithout Interaction term
linM = train %>%
    select(RNA_exp_log, UP_element, Minus35, Spacer, Minus10, Background) %>%
    lm(formula=RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=., model = T) 

test$predicted_exp_lm <- predict(linM, test) #predict expression of test data

corr <- summary(lm(RNA_exp_log ~ predicted_exp_lm, test))$r.square

pal <- wes_palette("Zissou", 8, type = "continuous")
test$Minus35  <- with(test, reorder(Minus35,RNA_exp_12, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_12, median)))

ggplot(test, aes(10^predicted_exp_lm, 10^RNA_exp_log, color = Minus35)) +
    geom_point(alpha = .4, aes(color = Minus35)) +
    annotate("text", x = .05, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
    scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
    scale_x_log10(limits = c(0.01,30)) + scale_y_log10(limits = c(0.01,30)) + 
    annotation_logticks() + 
    xlab('Predicted Expression') + ylab('Expression') +
    ggtitle('log-linear no interaction term') 
ggsave('log-linear_noint.pdf')

# Supp4D, neural network with only 5% of training data
smp_size <- floor(0.05 * nrow(split_exp))
set.seed(124)
train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)

train <- split_exp[train_ind, ]
test <- split_exp[-train_ind, ]

set.seed(124)
fit <- nnet(RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=train, size=10, maxit=300, linout=T, decay=0.01)

# make predictions
test$predicted_exp_NN <- predict(fit, test, type="raw")

corr <- summary(lm(RNA_exp_log ~ predicted_exp_NN, test))$r.square

pal <- wes_palette("Zissou", 8, type = "continuous")
test$Minus35  <- with(test, reorder(Minus35,RNA_exp_1, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_1, median)))

ggplot(test, aes(10^predicted_exp_NN, 10^RNA_exp_log, color = Minus35)) + 
    geom_point(alpha = .4) +
    annotate("text", x =0.1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
    scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) + 
    scale_x_log10(limits = c(0.01,30)) + scale_y_log10(limits = c(0.01,30)) + 
    annotation_logticks() + 
    xlab('Predicted Expression') +
    ylab('Expression') +
    ggtitle('Neural Network on .05')
ggsave('5percent_NN.pdf')

# Supplementary Figure 5, Neural Network schematic

install.packages('devtools')
library(devtools)
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
plot.nnet(fit)

# Supplementary Figure 6
# This figure is in two parts
pal <- wes_palette("Zissou",8, type = "continuous")

split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_1, median))
legend_ord <- levels(with(split_exp, reorder(Minus35,RNA_exp_1, -median)))

split_exp$Minus10 <- with(split_exp, reorder(Minus10,RNA_exp_1, median))
legend_ord <- levels(with(split_exp, reorder(Minus10,-RNA_exp_1, median)))

# 326xUP Log2 (TOP)
split_exp %>% 
    filter(UP_element == 'gourse-326fold-up') %>% 
    select(RNA_exp_12, Minus35, Minus10) %>% filter(RNA_exp_12 > 0) %>%
    group_by_(.dots=c("Minus35", "Minus10")) %>% 
    summarize(RNA_exp_12 = median(RNA_exp_12)) %>% 
    ungroup() %>%
    ggplot(., aes(x=Minus10, y= Minus35)) + 
    geom_tile(aes(fill = log2(RNA_exp_12))) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradientn(colors = pal, name = 'Log2(Expression)', limits = c(-4,3.3)) +
    labs(title = 'Combinatorial Expression of Core Promoter Elements (326xUP Element)', 
         x = '-10 Variant', y = '-35 Variant')
ggsave("Combinatorial_326xUP_log2.pdf")

# NO UP (Bottom)
split_exp %>% 
    filter(UP_element == 'noUP') %>% 
    select(RNA_exp_12, Minus35, Minus10) %>% 
    filter(RNA_exp_12 > 0) %>%
    group_by_(.dots=c("Minus35", "Minus10")) %>% 
    summarize(RNA_exp_12 = median(RNA_exp_12)) %>% 
    ungroup() %>%
    ggplot(., aes(x=Minus10, y= Minus35)) + 
    geom_tile(aes(fill = log2(RNA_exp_12))) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradientn(colors = pal, name = 'Log2(Expression)', limits = c(-4,3.3)) +
    labs(title = 'Combinatorial Expression of Core Promoter Elements (326xUP Element)', 
         x = '-10 Variant', y = '-35 Variant')
ggsave("Combinatorial_noUP_log2.pdf")

# Supplementary Figure 7 
pal <- wes_palette("Zissou",8, type = "continuous")

split_exp$Minus10 <- with(split_exp, reorder(Minus10,RNA_exp_12, median))
legend_ord <- levels(with(split_exp, reorder(Minus10,-RNA_exp_12, median)))

a <- split_exp %>% 
    filter(UP_element == 'gourse-326fold-up') %>% 
    select(RNA_exp_12, Minus35, Minus10, Spacer, Background) %>%
    filter(RNA_exp_12 > -0) %>%
    ungroup()

b <-split_exp %>% 
    filter(UP_element == 'noUP') %>% 
    select(RNA_exp_12, Minus35, Minus10, Spacer, Background) %>%
    filter(RNA_exp_12 > -0) %>%
    ungroup()

inner_join(a,b, by = c('Minus10', 'Minus35', 'Background', 'Spacer')) %>%
    transform(Fold_UP = RNA_exp_12.x/RNA_exp_12.y) %>% 
    ggplot(., aes(x=RNA_exp_12.y, y=Fold_UP)) + 
    geom_point(aes(color = Minus10)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    annotation_logticks(sides = 'l') +
    scale_color_manual(values = pal, name = '-10 Variant', breaks = legend_ord) +
    ylim(0,25) + scale_y_log10() + geom_hline(yintercept=1, color = 'black') +
    labs(title = 'Increase Due to UP-Element', x = 'Base Promoter Strength', 
         y = 'Increase in Expression Due to UP Element')
ggsave("UPIncrease_log.pdf")


# Supplementary Figure 8
split_exp <- mutate(split_exp, RNA_exp_log = log10(RNA_exp_12)) #Log transformed data for lm

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

pal <- wes_palette("Zissou", 8, type = "continuous")
test$Background  <- with(test, reorder(Background,RNA_exp_12, median))
legend_ord <- levels(with(test, reorder(Background,-RNA_exp_12, median)))

ggplot(test, aes(10^predicted_exp_NN, 10^RNA_exp_log, color = Background)) + 
    geom_point(alpha = .4) +
    annotate("text", x =0.1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T, size = 6) +
    scale_color_manual(values = pal, name = 'Background', breaks = legend_ord) + 
    scale_x_log10(limits = c(0.01,30)) + scale_y_log10(limits = c(0.01,30)) + 
    annotation_logticks() + 
    xlab('Predicted Expression') +
    ylab('Expression') +
    ggtitle('Using Neural Network to Predict Variant Expression')
ggsave('SupplementaryFigure6.png')


# Supplementary Figure 9
pal <- wes_palette("Zissou",8, type = "continuous")

split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_12, median))
legend_ord <- levels(with(split_exp, reorder(Minus35,RNA_exp_12, -median)))

split_exp$Minus10 <- with(split_exp, reorder(Minus10,RNA_exp_12, median))
legend_ord <- levels(with(split_exp, reorder(Minus10,-RNA_exp_12, median)))

# ECK726, TOP
split_exp %>% 
    filter(Spacer == 'ECK125137726') %>% 
    select(RNA_exp_12, Minus10, Minus35) %>% 
    filter(RNA_exp_12 > 0) %>%
    group_by_(.dots=c("Minus10", "Minus35")) %>% 
    summarize(RNA_exp_12 = median(RNA_exp_12)) %>% 
    ungroup() %>% 
    ggplot(., aes(x=Minus10, y=Minus35)) + 
    geom_tile(aes(fill = log10(RNA_exp_12))) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradientn(colors = pal, name = 'Log10(Expression)', limits = c(-1.2,1)) +
    labs(title = 'ECK726', x = '-10 Variant', y = '-35 Variant')
ggsave('Sup7_ECK726.pdf')

# ECK405, Bottom
split_exp %>% 
    filter(Spacer == 'ECK125137405')%>% 
    select(RNA_exp_12, Minus10, Minus35) %>% 
    filter(RNA_exp_12 > 0) %>%
    group_by_(.dots=c("Minus10", "Minus35")) %>% 
    summarize(RNA_exp_12 = median(RNA_exp_12)) %>% 
    ungroup() %>%
    ggplot(., aes(x=Minus10, y=Minus35)) + 
    geom_tile(aes(fill = log10(RNA_exp_12))) +  
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradientn(colors = pal, name = 'Log10(Expression)', limits = c(-1.2,1)) +
    labs(title = 'ECK405', x = '-10 Variant', y = '-35 Variant')
ggsave('Sup7_ECK405.pdf')
