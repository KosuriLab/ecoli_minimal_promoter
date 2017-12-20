library("wesanderson")
library("ggplot2")
library("dplyr")

pal <- wes_palette("Zissou", 8, type = "continuous")

split_exp <- read.table('../processed_data/rlp5Min_SplitVariants.txt', header = T)

split_exp$Minus35  <- with(split_exp, reorder(Minus35,RNA_exp_average, median))
split_exp$Minus10 <- with(split_exp, reorder(Minus10,-RNA_exp_average, median))
split_exp$Background <- with(split_exp, reorder(Background,RNA_exp_average, median))
split_exp$Spacer <- with(split_exp, reorder(Spacer,RNA_exp_average, median))

#Figure 3A
split_exp$UP_element <- with(split_exp, reorder(UP_element,RNA_exp_average, max))

ggplot(filter(split_exp, Minus10 == 'consensus10', 
              Minus35 == 'consensus35', 
              Background != 'bg4323949:4324099', 
              Background != 'bg977040:977190'), aes(x = Spacer, y = Background)) + 
  geom_tile(aes(fill = log10(RNA_exp_average))) + 
  scale_fill_gradientn(colors = pal, limits = c(-1.4,1.6)) + 
  facet_wrap(~UP_element) +
  theme(panel.background = element_rect(fill = "gray80"), 
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 6), 
        axis.text.x = element_text(size = 6, angle = -90), 
        strip.text = element_text(size = 6),
        panel.spacing = unit(0.40, "lines"))
ggsave('../figs/Figure3A.pdf')

#Figure 3B
ggplot(filter(split_exp, 
              UP_element == 'gourse-136fold-up', 
              Background != 'bg4323949:4324099', 
              Background != 'bg977040:977190'), aes(x = Spacer, y = Background)) + 
  geom_tile(aes(fill = log10(RNA_exp_average))) + 
  scale_fill_gradientn(colors = pal, limits = c(-1.4,1.6)) + 
  facet_grid(Minus10 ~ Minus35) +
  theme(panel.background = element_rect(fill = "gray80"), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank(), 
        strip.text = element_blank(),
        panel.spacing = unit(0.10, "lines"))
ggsave('../figs/Figure3B.pdf')


                                                                                                