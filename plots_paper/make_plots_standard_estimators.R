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

experiment.df <- read.csv("../data/standard_estimators.csv")

estimator.rename <- function(estimator) {
  if (estimator == "lin.out.intercept") {
    return("OLS Intercept")
  } else if (estimator == "lin.out.cf1.g.est") {
    return("1-fold G")
  } else if (estimator == "lin.out.cf2.g.est") {
    return("2-fold G")
  } else if (estimator == "nlin.out.cf1.g.est") {
    return("1-fold G (misspecified)")
  } else if (estimator == "nlin.out.cf2.g.est") {
    return("2-fold G (misspecified)")
  } else if (estimator == "lin.out.cf1.g.reg.est") {
    return("1-fold G (ridge)")
  } else if (estimator == "lin.out.cf2.g.reg.est") {
    return("2-fold G (ridge)")
  } else if (estimator == "aipw.est_prop.lin_out.cf1.est") {
    return("1-fold AIPW")
  } else if (estimator == "aipw.est_prop.lin_out.cf2.est") {
    return("2-fold AIPW")
  } else if (estimator == "aipw.est_prop.lin_out.cf3.est") {
    return("3-fold AIPW")
  } else if (estimator == "ipw.est_prop.lin_out.cf1.est") {
    return("1-fold IPW")
  } else if (estimator == "ipw.est_prop.lin_out.cf2.est") {
    return("2-fold IPW")
  } else if (estimator == "lin.out.intercept") {
    return("OLS Intercept")
  } else if (estimator == "aipw.est_prop.nlin_out.cf1.est") {
    return("1-fold AIPW (misspecified)")
  } else if (estimator == "aipw.est_prop.nlin_out.cf2.est") {
    return("2-fold AIPW (misspecified)")
  } else if (estimator == "aipw.est_prop.nlin_out.cf3.est") {
    return("3-fold AIPW (misspecified)")
  } else if (estimator == "aipw.est_prop.reg_lin_out.cf1.est") {
    return("1-fold AIPW (ridge)")
  } else if (estimator == "aipw.est_prop.reg_lin_out.cf2.est") {
    return("2-fold AIPW (ridge)")
  } else if (estimator == "aipw.est_prop.reg_lin_out.cf3.est") {
    return("3-fold AIPW (ridge)")
  } else {
    return(estimator)
  }
}
estimator.rename = Vectorize(estimator.rename)
colnames(experiment.df)

plotting.df <-
  select(experiment.df,!c(seed,p,X)) %>%
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

### Well Specified, comparison of Outcome modeling, AIPW, IPW ###
# Define the desired order of categories
category_order <- c("1-fold G",
                    "2-fold G",
                    "2-fold AIPW",
                    "3-fold AIPW",
                    "2-fold IPW")

# Plots specs
# # Uncomment for presentation plots
# my_theme <- theme(legend.position = "bottom",
#                   legend.title = element_blank(),
#                   axis.text = element_text(size = 13),
#                   axis.title = element_text(size = 17),
#                   legend.text = element_text(size = 13))
# linewidth <- 1
# point_size <- 2
# errbar_width <- .05
# errbar_linewidth <- 1

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

# Create the plot
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
  geom_hline(yintercept=0,linetype="dashed",linewidth=linewidth) +
  scale_color_manual(values = scales::hue_pal()(length(category_order)),  # Use default color palette
                     breaks = category_order,
                     labels = category_order) +
  scale_shape_manual(values = 1:length(category_order),  # Set shape values
                     breaks = category_order,
                     labels = category_order) +
  scale_linetype_manual(values = 1:length(category_order),  # Set linetype values
                        breaks = category_order,
                        labels = category_order) +
  scale_x_log10() +
  theme_bw() +
  labs(y = "Mean of estimate") +
  my_theme +
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         shape = guide_legend(nrow = 1, byrow = TRUE),
         linetype = guide_legend(nrow = 1, byrow = TRUE),
         override.aes = list(shape = 16))
mean.plot

# var plot
var.plot <- ggplot(data=filter(plotting.df,
                               estimator %in% category_order,
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
         linetype = guide_legend(nrow = 2, byrow = TRUE),
         override.aes = list(shape = 16))
var.plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
my_legend <- g_legend(mean.plot)

pdf("../fig_paper/well_specified_OLS_outcome_comparison.pdf",width=8,height=5)
grid.arrange(arrangeGrob(mean.plot + theme(legend.position="none"),
                         var.plot + theme(legend.position="none"),
                         nrow=1), my_legend,
             nrow = 2, heights=c(10, 2))
dev.off()

pdf("../fig_paper/well_specified_OLS_outcome_comparison_MEAN.pdf",width=4,height=4)
mean.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/well_specified_OLS_outcome_comparison_VAR.pdf",width=4,height=4)
var.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/well_specified_OLS_outcome_comparison_LEGEND.pdf",width=6,height=.5)
grid.draw(my_legend)
dev.off()


### Mis-specified, comparison of Outcome modeling and AIPW ###
# Define the desired order of categories
category_order <- c("1-fold G (misspecified)",
                    "2-fold G (misspecified)",
                    "2-fold AIPW (misspecified)",
                    "3-fold AIPW (misspecified)")

# Create the plot
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
  geom_hline(yintercept=0,linetype="dashed",linewidth=linewidth) +
  scale_color_manual(values = scales::hue_pal()(length(category_order)),  # Use default color palette
                     breaks = category_order,
                     labels = category_order) +
  scale_shape_manual(values = 1:length(category_order),  # Set shape values
                     breaks = category_order,
                     labels = category_order) +
  scale_linetype_manual(values = 1:length(category_order),  # Set linetype values
                        breaks = category_order,
                        labels = category_order) +
  scale_x_log10() +
  theme_bw() +
  labs(y = "Mean of estimate") +
  my_theme +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE),
         override.aes = list(shape = 16))
mean.plot

# var plot
var.plot <- ggplot(data=filter(plotting.df,
                               estimator %in% category_order,
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
         linetype = guide_legend(nrow = 2, byrow = TRUE),
         override.aes = list(shape = 16))
var.plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
my_legend <- g_legend(mean.plot)

pdf("../fig_paper/misspecified_OLS_outcome_comparison.pdf",width=8,height=5)
grid.arrange(arrangeGrob(mean.plot + theme(legend.position="none"),
                         var.plot + theme(legend.position="none"),
                         nrow=1), my_legend,
             nrow = 2, heights=c(10, 2))
dev.off()

pdf("../fig_paper/misspecified_OLS_outcome_comparison_MEAN.pdf",width=4,height=4)
mean.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/misspecified_OLS_outcome_comparison_VAR.pdf",width=4,height=4)
var.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/misspecified_OLS_outcome_comparison_LEGEND.pdf",width=6,height=.5)
grid.draw(my_legend)
dev.off()


### Mis-specified, comparison of Outcome modeling and AIPW ###
# Define the desired order of categories
category_order <- c("1-fold G (ridge)",
                    "2-fold G (ridge)",
                    "2-fold AIPW (ridge)",
                    "3-fold AIPW (ridge)")

# Create the plot
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
  geom_hline(yintercept=0,linetype="dashed",linewidth=linewidth) +
  scale_color_manual(values = scales::hue_pal()(length(category_order)),  # Use default color palette
                     breaks = category_order,
                     labels = category_order) +
  scale_shape_manual(values = 1:length(category_order),  # Set shape values
                     breaks = category_order,
                     labels = category_order) +
  scale_linetype_manual(values = 1:length(category_order),  # Set linetype values
                        breaks = category_order,
                        labels = category_order) +
  scale_x_log10() +
  theme_bw() +
  labs(y = "Mean of estimate") +
  my_theme +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE),
         override.aes = list(shape = 16))
mean.plot

# var plot
var.plot <- ggplot(data=filter(plotting.df,
                               estimator %in% category_order,
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
         linetype = guide_legend(nrow = 2, byrow = TRUE),
         override.aes = list(shape = 16))
var.plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
my_legend <- g_legend(mean.plot)

pdf("../fig_paper/well_specified_ridge_outcome_comparison.pdf",width=8,height=5)
grid.arrange(arrangeGrob(mean.plot + theme(legend.position="none"),
                         var.plot + theme(legend.position="none"),
                         nrow=1), my_legend,
             nrow = 2, heights=c(10, 2))
dev.off()

pdf("../fig_paper/well_specified_ridge_outcome_comparison_MEAN.pdf",width=4,height=4)
mean.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/well_specified_ridge_outcome_comparison_VAR.pdf",width=4,height=4)
var.plot + theme(legend.position="none")
dev.off()

pdf("../fig_paper/well_specified_ridge_outcome_comparison_LEGEND.pdf",width=6,height=.5)
grid.draw(my_legend)
dev.off()
