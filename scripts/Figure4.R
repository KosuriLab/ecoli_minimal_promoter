#install.packages('nnet')
#install.packages('dplyr')
#install.packages('ggplot2')
#install.packages('wesanderson')
#install.packages('nnet')
#install.packages('cowplot')
#install.packages('reshape2')

library("ggplot2")
library("dplyr")
library("wesanderson")
library("reshape2")
require(cowplot)
library(nnet)

options(stringsAsFactors = F)

split_exp <- read.table('rlp5Min_SplitVariants.txt', header = T)
split_exp <- mutate(split_exp, RNA_exp_log = log10(RNA_exp_average))

#Sample training and test data

smp_size <- floor(0.5 * nrow(split_exp))

set.seed(123)
train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)

train <- split_exp[train_ind, ]
test <- split_exp[-train_ind, ]



#Figure 4A, Predict expression using linear model trained on 50% of log-transformed data


#train model to predict expression based on each individual element and combinations of -10 and -35 elements
linM = train %>%
  select(RNA_exp_log, UP_element, Minus35, Spacer, Minus10, Background) %>%
  lm(formula=RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background + Minus10:Minus35, data=., model = T) 

test$predicted_exp_lm <- predict(linM, test) #predict expression of test data

corr <- summary(lm(RNA_exp_log ~ predicted_exp_lm, test))$r.square



pal <- wes_palette("Zissou", 8, type = "continuous")
test$Minus35  <- with(test, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_average, median)))


ggplot(test, aes(10^predicted_exp_lm, 10^RNA_exp_log, color = Minus35)) +
  geom_point(alpha = .4, aes(color = Minus35)) +
  annotate("text", x = .1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
  scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) +
  scale_x_log10(limits = c(0.05,30)) + scale_y_log10(limits = c(0.05,30)) + annotation_logticks() + 
  xlab('Predicted Expression') + ylab('Expression') +
  ggtitle('Using Log-Linear Model to Predict Variant Expression') 
ggsave('Figure4A.pdf')



#Figure 4B

#We use ANOVA to determine the relative contribution of each element to expression

anova(linM)

#Explained variance as a percentage of the sum of squared deviation

#FIRST TIME
pie_linM <- data.frame(
  Element = c("Minus10", "Minus35", 'UP element', 'Spacer', 'Background', 'Minus10:Minus35', 'Residual' ),
  value = c(413.14, 648.09,  47.94,  74.89, 21.82, 380.18,368.84) #Mean squares for 
)


pie_linM$Element  <- with(pie_linM, reorder(Element,value, median))
pie_linM <- mutate(pie_linM, temp = value/sum(value))

ggplot(pie_linM, aes(x=Element, y=temp))+ coord_flip() +
  geom_bar(width = 1, stat = "identity", fill='#3B9AB2', color = 'Black') + 
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Explained Variance by Element", y = 'Explained Variance', x='Element')
ggsave('Figure4B.pdf')


#Figure 4C

#train neural network on all sequence elements
set.seed(123)
fit <- nnet(RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=train, size=10, maxit=300, linout=T, decay=0.01)

# make predictions
test$predicted_exp_NN <- predict(fit, test, type="raw")

corr <- summary(lm(RNA_exp_log ~ predicted_exp_NN, test))$r.square

pal <- wes_palette("Zissou", 8, type = "continuous")
test$Minus35  <- with(test, reorder(Minus35,RNA_exp_average, median))
legend_ord <- levels(with(test, reorder(Minus35,-RNA_exp_average, median)))

ggplot(test, aes(10^predicted_exp_NN, 10^RNA_exp_log, color = Minus35)) + 
  geom_point(alpha = .4) +
  annotate("text", x =0.1, y = 10, label = paste('R^2==', signif(corr, 3)), parse = T) +
  scale_color_manual(values = pal, name = '-35 Variant', breaks = legend_ord) + 
  scale_x_log10(limits = c(0.05,30)) + scale_y_log10(limits = c(0.05,30)) + annotation_logticks() + 
  xlab('Predicted Expression') +
  ylab('Expression') +
  ggtitle('Using Neural Network to Predict Variant Expression')
ggsave('Figure4C.pdf')


#Figure 4D, cross-validate  Get R^2 for various subsets
a= data.frame()  
corrz <- list() # new empty list

set.seed(123)
for (i in c(.05,.1, .2,.3,.4, .5, .6, .7, .8, .9)) {
  smp_size <- floor(i * nrow(split_exp))
  
  train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)
  train <- split_exp[train_ind, ]
  test <- split_exp[-train_ind, ]
  
  fit <- nnet(RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=train, size=10, maxit=300, linout=T, decay=0.01)
  test$predicted_exp_NN <- predict(fit, test, type="raw")
  
  p1 = summary(lm(RNA_exp_log ~ predicted_exp_NN, test))$r.square
  corrz[[as.character(i)]] <- p1
}

a<- as.numeric(corrz) #DO this AFTER THE FIRST LOOP OR YOU'LL GET an error

set.seed(123)
for (n in c(1:9)) {
for (i in c(.05,.1, .2,.3,.4, .5, .6, .7, .8, .9)) {
  smp_size <- floor(i * nrow(split_exp))

  train_ind <- sample(seq_len(nrow(split_exp)), size = smp_size)
  train <- split_exp[train_ind, ]
  test <- split_exp[-train_ind, ]
  
  fit <- nnet(RNA_exp_log ~ UP_element + Minus35 + Spacer + Minus10 + Background, data=train, size=10, maxit=300, linout=T, decay=0.01)
  test$predicted_exp_NN <- predict(fit, test, type="raw")

  p1 = summary(lm(RNA_exp_log ~ predicted_exp_NN, test))$r.square
  corrz[[as.character(i)]] <- p1
}

  a <- cbind(a, sapply(corrz, function(x){as.numeric(x[1])}))
}

regression <- as.data.frame(a)
regression['subset']<- row.names(regression)
melt<- melt(regression, by = 'subset')

ggplot(melt, aes(x=subset, y=value)) + stat_boxplot(geom = "errorbar", width = 0.25) + 
  geom_boxplot(lwd =.5, fill = '#3B9AB2') +
  geom_jitter() + 
  ylim(.75,1) +
  labs(x='Proportion of Data Trained On', y='Correlation between Predicted and Actual Expression', title = '10-Fold Cross Validation with Variable Proportion of Training Data')
ggsave('Figure4D.pdf')
