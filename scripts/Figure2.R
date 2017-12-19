# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("cowplot")
# install.packages("wesanderson")
# install.packages("tidyr")
# install.packages("reshape2")


library("ggplot2")
library("dplyr")
library("wesanderson")
library("reshape2")
library('tidyr')
require(cowplot)


min_MOPS <- read.table("./revLP5_Min_MOPS_glu_expression.txt", header = T)

# Take promoter variant names and assign columns indicating which element variants it has
variant_stats <- read.table("./variant_statistics.txt", header = T, sep = "\t")

# First, rename UP elements so we can split by _ 
variant_stats$name <- gsub("gourse_136fold_up", "gourse-136fold-up", variant_stats$name)
variant_stats$name <- gsub("gourse_326fold_up", "gourse-326fold-up", variant_stats$name)
min_MOPS$name <- gsub("gourse_136fold_up", "gourse-136fold-up", min_MOPS$name)
min_MOPS$name <- gsub("gourse_326fold_up", "gourse-326fold-up", min_MOPS$name)

# Now split variant names 
split_exp <- variant_stats %>% 
    separate(col = 'name', into = c("UP_element", "Minus35", "Spacer", "Minus10", 'Background'),
             sep = "_", remove = F) %>%
    left_join(., min_MOPS, by = 'name') 
# The 'too few values' warning is referring to the negative and positive controls 
# within our library, which do not have these elements and so get removed

split_exp <- split_exp[!is.na(split_exp$Background),]
split_exp <- split_exp[!is.na(split_exp$RNA_exp_1),]

write.table(split_exp, "rlp5Min_SplitVariants.txt", quote = F, row.names = F)

# Figure 2C, distribution of the number of mapped and integrated barcodes

ggplot(min_MOPS, aes(num_barcodes)) + 
    geom_density(alpha=.2, fill = "#3B9AB2", color = "#3B9AB2") +
    geom_density(aes(num_barcodes_integrated), alpha=.2, fill = "#F21A00", color = "#F21A00") +
    scale_x_log10() + annotation_logticks(sides = 'b') +
    annotate("text", x=7, y=1.35, label = paste(' Integrated Median =', median(min_MOPS$num_barcodes_integrated), 'barcodes')) +
    annotate("text", x=7, y=1.43, label = paste('Mapped Median =', median(min_MOPS$num_barcodes), 'barcodes')) +
    labs(x='Number of Mapped Barcodes', y='Density', title = 'Distribution of the Number of Barcodes Per Variant')

ggsave('Figure2C.pdf')

# Figure 2E
# Plot expression correlations with negative controls and consensus promoters highlighted
test <- subset(min_MOPS, grepl("neg_control", min_MOPS$name))
corr <- summary(lm(RNA_exp_1 ~ RNA_exp_2, min_MOPS))$r.squared
summary(lm(RNA_exp_1 ~ RNA_exp_2, min_MOPS)) # p < 2.2e-16

ggplot(min_MOPS, aes(RNA_exp_1, RNA_exp_2)) + geom_point(alpha = .4) +
    annotate("text", x =.05, y = 5, label = paste('R^2==', signif(corr, 3)), parse = T) + 
    annotation_logticks() + scale_x_log10(breaks = c(1,10,100)) + scale_y_log10(breaks = c(1,10,100)) + 
    xlab('MOPS + .2% Glu Expression Rep1') + ylab('MOPS + .2% Glu Expression Rep2') +
    ggtitle('Comparing Promoter Expression Between Technical Replicates') +
    geom_point(data=test, aes(RNA_exp_1, RNA_exp_2), color = "firebrick1") + 
    geom_point(data=filter(split_exp, Minus10 == 'consensus10' & Minus35 == 'consensus35'),
               aes(RNA_exp_1,RNA_exp_2), color = 'deepskyblue')
ggsave('Figure2E.pdf')


