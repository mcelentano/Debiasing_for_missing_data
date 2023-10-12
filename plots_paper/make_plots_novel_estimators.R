# run from the plots_paper directory

library(tidyverse)
library(dplyr)
library(reshape)
library(viridis)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(viridis)

# Load and wrangle data
all_data.df <- read.csv("../data/novel_estimators.csv")

estimator.rename <- function(estimator) {
  if (estimator == "mu.y.hat") {
    return("Ridge")
  } else if (estimator == "mu.y.hat.ipw") {
    return("Ridge (IPW)")
  } else if (estimator == "mu.y.hat.d.naive") {
    return("Debiased ridge (naive)")
  } else if (estimator == "mu.y.hat.d.ipw") {
    return("Debiased ridge (IPW)")
  } else if (estimator == "mu.y.hat.d.cfd") {
    return("Oracle ASCW")
  } else if (estimator == "mu.y.hat.d.emp.moment") {
    return("Empirical SCA")
  } else {return(estimator)}
}
estimator.rename = Vectorize(estimator.rename)


plotting.df <-
  select(all_data.df,!c(seed,p,sigma,theta.d.0,theta.y.0,lambda.y,lambda.d,mu.y,mu.d,gamma.y,gamma.d,pi.bar,alpha.1,alpha.2)) %>%
  melt(id=c("n")) %>%
  group_by(n,variable) %>%
  summarize(
    mean = mean(value,na.rm=TRUE), # why do I get some na's for small samples?
    var = var(value,na.rm=TRUE),
    mean_se = sd(value,na.rm=TRUE)/sqrt(n()),
    var_se = sqrt( var( (value - mean(value))^2 ) / (n()-2) )
  ) %>%
  ungroup %>% as.data.frame %>%
  dplyr::rename(estimator=variable) %>%
  melt(id=c("n","estimator","mean_se","var_se")) %>%
  dplyr::rename(metric=variable,mean=value) %>%
  mutate(se = if_else(metric=="mean",mean_se,var_se)) %>% 
  mutate(estimator = as.character(estimator),
         estimator = estimator.rename(estimator)) %>% 
  select(c(n,estimator,metric,mean,se))
plotting.df$estimator

# mean plot
mean.plot <- ggplot(data = filter(plotting.df,
                                  estimator %in% c("Ridge",
                                                   "Ridge (IPW)",
                                                   "Debiased ridge (naive)",
                                                   "Debiased ridge (IPW)",
                                                   "Oracle ASCW",
                                                   "Empirical SCA"),
                                  metric == "mean")) +
  geom_line(mapping = aes(x = n, y = mean,
                          color = estimator),
            linewidth = 1) +
  geom_point(mapping = aes(x = n, y = mean,
                           color = estimator,
                           shape = estimator),
             size = 2) +
  geom_errorbar(aes(x = n,
                    ymin = mean - 1.96 * se, ymax = mean + 1.96 * se,
                    color = estimator),
                width = .05, linewidth = 1) +
  scale_x_log10() +
  theme_bw() +
  labs(y = "Mean of estimate") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
mean.plot

# Define the desired order of categories
category_order <- c("Ridge", "Debiased ridge (naive)",
                    "Debiased ridge (IPW)", "Oracle ASCW",
                    "Empirical SCA")

# Uncomment for paper plots
my_theme <- theme(legend.position = "bottom",
                  legend.title = element_blank(),
                  axis.text = element_text(size = 10),
                  axis.title = element_text(size = 13),
                  legend.text = element_text(size = 10))
linewidth <- .5
point_size <- 2
errbar_width <- .05
errbar_linewidth <- .5

mean.plot <- ggplot(data = filter(plotting.df,
                                  estimator %in% category_order,
                                  metric == "mean")) +
  geom_line(mapping = aes(x = n, y = mean,
                          color = estimator,
                          linetype = estimator),
            linewidth = linewidth) +  # Remove alpha setting
  geom_point(mapping = aes(x = n, y = mean,
                           color = estimator,
                           shape = estimator),
             size = point_size) +
  geom_errorbar(aes(x = n,
                    ymin = mean - 1.96 * se, ymax = mean + 1.96 * se,
                    color = estimator),
                width = errbar_width, linewidth = errbar_linewidth) +
  scale_color_manual(values = scales::hue_pal()(length(category_order)),  # Use default color palette
                     breaks = category_order,
                     labels = category_order) +
  scale_shape_manual(values = 1:length(category_order),  # Set shape values
                     breaks = category_order,
                     labels = category_order) +
  scale_linetype_manual(values = 1:length(category_order),  # Set linetype values
                        breaks = category_order,
                        labels = category_order) +
  geom_hline(yintercept=0,linetype="dashed",linewidth=1) +
  scale_x_log10() +
  theme_bw() +
  labs(y = "Mean of estimate") +
  my_theme +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE),
         override.aes = list(shape = 16))
mean.plot

# var plot
var.plot <- ggplot(data=filter(plotting.df,
                               estimator%in%category_order,
                               metric=="var")) +
  geom_line(mapping=aes(x=n,y=mean,
                        color=estimator,
                        linetype = estimator),
            linewidth=linewidth) +
  geom_point(mapping=aes(x=n,y=mean,
                         color=estimator,
                         shape=estimator),
             size=point_size) +
  geom_errorbar(
    aes(x=n,
        ymin=mean-1.96*se,ymax=mean+1.96*se,
        color=estimator),
    width=errbar_width,linewidth=errbar_linewidth
  ) +
  scale_color_manual(values = scales::hue_pal()(length(category_order)),  # Use default color palette
                     breaks = category_order,
                     labels = category_order) +
  scale_shape_manual(values = 1:length(category_order),  # Set shape values
                     breaks = category_order,
                     labels = category_order) +
  scale_linetype_manual(values = 1:length(category_order),  # Set linetype values
                        breaks = category_order,
                        labels = category_order) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  labs(y="Variance of estimate") +
  my_theme +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE),
         override.aes = list(shape = 16))
var.plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
my_legend <- g_legend(mean.plot)

pdf("../fig_paper/debiasing.pdf",width=8,height=5)
grid.arrange(arrangeGrob(mean.plot + theme(legend.position="none"),
                         var.plot + theme(legend.position="none"),
                         nrow=1), my_legend,
             nrow = 2, heights=c(10, 2))
dev.off()

pdf("../fig_paper/debiasing_MEAN.pdf",width=4,height=4)
mean.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/debiasing_VAR.pdf",width=4,height=4)
var.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/debiasing_LEGEND.pdf",width=6,height=.8)
grid.draw(my_legend)
dev.off()
