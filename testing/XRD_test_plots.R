library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(gtable)
library(egg)
library(viridis)
library(stringr) 

theme_set(
  theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
)

PE_4P <- read.table('4-6P_pattern_errors.txt', sep='', header=TRUE, row.names=NULL)
PE_8P <- read.table('8-12P_pattern_errors.txt', sep='', header=TRUE, row.names=NULL)
PE_4P$PC <- '4-6%'
PE_8P$PC <- '8-12%'
PE <- rbind(PE_4P, PE_8P)
PE$MAE <- PE$MAE/max(PE$MAE)
PE$MSE <- PE$MSE/max(PE$MSE)
PE$Clark <- PE$Clark/max(PE$Clark)
PE$PS <- PE$PS/max(PE$PS)
PE$niter <- PE$niter/max(PE$niter)

PE <- melt(PE, id.vars=c('name','EF', 'PC'), measure.vars=c('MAE', 'MSE', 'Clark', 'PS', 'niter'))

medians <- PE %>%
  group_by(EF, variable, PC) %>%
  summarise_at(vars(value), funs(median(., na.rm=TRUE)))

names(medians)[names(medians) == 'value'] <- 'median'

medians <- medians %>%
  group_by(variable, PC) %>%
  mutate(median_rank=order(order(median, decreasing=FALSE)))

medians[medians$median_rank > 3, ]$median_rank <- '> 3'
medians$median_rank <- factor(medians$median_rank)

PE <- merge(PE, medians, by=c('variable','EF','PC'))
PE <- subset(PE, variable != 'Clark')

p1 <- ggplot() +
  geom_boxplot(data=PE, aes(x=EF, y=value, fill=median_rank),
               outlier.shape=4, outlier.color='black') +
  scale_fill_manual(values=c('white', 'goldenrod1', 'indianred1', 'deepskyblue1')) +
  facet_wrap(~ PC + variable, ncol=4) +
  coord_cartesian(ylim=c(0,0.5)) +
  labs(x='Error Function', y='Normalized Error Metric') +
  theme(axis.text.x=element_text(angle=90, hjust=1.0, vjust=0.5),
        strip.background=element_blank(),
        strip.text=element_blank())

ggsave('pattern_errors.tiff', plot=p1, width=10, height=5.5)

UC_4P <- read.table('4-6P_uc_errors.txt', sep='', header=TRUE, row.names=NULL)
UC_8P <- read.table('8-12P_uc_errors.txt', sep='', header=TRUE, row.names=NULL)
UC_4P$PC <- '4-6% Initial Perturbation'
UC_8P$PC <- '8-12% Initial Perturbation'
UC <- rbind(UC_4P, UC_8P)
UC$PE <- 100 * abs(UC$matched - UC$correct)/UC$correct
UC <- subset(UC, param %in% c('a','b','c'))

UC <- UC %>%
  group_by(name, EF, PC) %>%
  summarise_at(vars(PE), funs(mean(., na.rm=TRUE)))

shapes=c(1,2,3,4,5,6,7,8,9,10)

p1 <- ggplot(UC, aes(x=EF, y=PE)) +
  geom_boxplot(alpha=0.25, outlier.shape=NA) + 
  geom_point(aes(shape=name, color=name)) +
  scale_shape_manual(values=shapes) +
  facet_wrap(~PC) +
  labs(x='Error Function', y='Mean Lattice Parameter Percent Error', shape='Material', color='Material') +
  theme(axis.text.x=element_text(angle=90, hjust=1.0, vjust=0.5))

ggsave('uc_errors.tiff', plot=p1, width=7.0, height=4.0)

quantile(subset(UC, EF=='city_block' & PC == '4-6% Initial Perturbation')$PE)
quantile(subset(UC, EF=='Euclidean' & PC == '4-6% Initial Perturbation')$PE)
quantile(subset(UC, EF=='squared_Euclidean' & PC == '4-6% Initial Perturbation')$PE)
quantile(subset(UC, EF=='MAE' & PC == '4-6% Initial Perturbation')$PE)
quantile(subset(UC, EF=='MSE' & PC == '4-6% Initial Perturbation')$PE)
quantile(subset(UC, EF=='probabilistic_symmetric' & PC == '4-6% Initial Perturbation')$PE)

quantile(subset(UC, EF=='city_block' & PC == '8-12% Initial Perturbation')$PE)
quantile(subset(UC, EF=='Euclidean' & PC == '8-12% Initial Perturbation')$PE)
quantile(subset(UC, EF=='squared_Euclidean' & PC == '8-12% Initial Perturbation')$PE)
quantile(subset(UC, EF=='MAE' & PC == '8-12% Initial Perturbation')$PE)




