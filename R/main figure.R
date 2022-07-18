# Title     : TODO
# Objective : TODO
# Created by: howanchung
# Created on: 2021/6/12

library(boot)
library(data.table)
library(magrittr)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(tibble)
library(gridExtra)
library(grid)
library(EpiEstim)
library(readxl)
library(viridis)
library(GGally)
library(network)
library(forcats)
library(dplyr)
library(arrow)
library(Cairo)
loadfonts()

TIME_CUT <- '2020-01-23'
FONTS <- 'STSong'
AGE_LABEL <- c('0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '>=60')

source('R/utils.R')
source('R/facet_wrap_custom.R')
source("R/data preparation.R")
source("R/statistical test.R")
test_results

####################COVID19 observed network visualization#########################
out_degree <- sna::degree(net, cmode = 'freeman')
col <- c('#ffeda0', viridis_pal(begin = 1, end = 0)(20)[c(1, 5, 7, 10, 15, 20)])
names(col) <- as.character(AGE_LABEL)
y <- plyr::mapvalues(network::get.vertex.attribute(net, 'age_interval'), from = 1:7,
                     to = AGE_LABEL)

set.seed(2021)
g1 <- ggnet2(net, color = y, palette = col, color.legend = "Age group",
             size = 'degree', size.min = 1, arrow.gap = 0.005, arrow.size = 3, edge.color = "household",
             edge.size = 0.6,
             legend.size = 30) +
    guides(size = FALSE,
           color = guide_legend(size = 1, override.aes = list(size = 5))) +
    #ggtitle("Overall observed transmission network") +
    theme(text = element_text(family = FONTS),
          legend.position = 'left',
          legend.key.size = unit(1, "cm"))


in_degree <- sna::degree(net, cmode = 'indegree')
out_degree <- sna::degree(net, cmode = 'outdegree')
nonsingletons <- which(in_degree != 0 | out_degree != 0)
out_degree <- table(out_degree[nonsingletons]) %>%
    as.data.table()

setnames(out_degree, c("outdegree", "count"))
out_degree[, outdegree := as.numeric(as.character(outdegree))]


ggplot(out_degree, aes(x = outdegree, y = count)) +
    geom_bar(stat = 'identity', width = 0.9) +
    geom_text(aes(label = count, y = count), vjust = -0.5, size = 5,
              family = FONTS) +
    scale_y_continuous(breaks = seq(0, 400, 100),
                       labels = seq(0, 400, 100),
                       limits = c(0, 410)) +
    labs(x = 'Out-degree', y = 'Frequency') +
    theme_bw() +
    annotation_custom(ggplotGrob(g1),
                      xmin = 2.5, xmax = 17, ymax = 405, ymin = 50) +
    theme(text = element_text(size = 30, family = FONTS),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

ggsave(filename = 'figure/Figure2.svg', width = 16, height = 12)
ggsave(filename = 'figure/Figure2.pdf', width = 16, height = 12, device=cairo_pdf)


#######
in_degree <- sna::degree(net, cmode = 'indegree')
out_degree <- sna::degree(net, cmode = 'outdegree')
nonsingletons <- which(in_degree != 0 | out_degree != 0)
mean(out_degree[nonsingletons])


degree_period <- data.table(id = 1:1349, out_degree = out_degree)
degree_period[, period := 2]
degree_period[id %in% (1:1349)[COVID19$origin_time <= TIME_CUT], period := 1]
degree_period <- degree_period[id %in% nonsingletons]

degree_period %>%
    group_by(period) %>%
    summarise(mean_outdegree = mean(out_degree))


########
set.seed(2021)
degree_ci1 <- boot.ci(boot(degree_period[period == 1, out_degree], bootMean, R = 10000), type = "perc")
set.seed(2021)
degree_ci2 <- boot.ci(boot(degree_period[period == 2, out_degree], bootMean, R = 10000), type = "perc")

degree_error <- data.table(group = c("I", "II"),
                           value = c(degree_ci1$t0, degree_ci2$t0),
                           low = c(degree_ci1$percent[4], degree_ci2$percent[4]),
                           upper = c(degree_ci1$percent[5], degree_ci2$percent[5]))
degree_error
########
set.seed(2021)
degree_ci1 <- boot.ci(boot(dist_period[period == 1, dist], bootMean, R = 10000), type = "perc")
set.seed(2021)
degree_ci2 <- boot.ci(boot(dist_period[period == 2, dist], bootMean, R = 10000), type = "perc")

dist_error <- data.table(group = c("I", "II"),
                         value = c(degree_ci1$t0, degree_ci2$t0),
                         low = c(degree_ci1$percent[4], degree_ci2$percent[4]),
                         upper = c(degree_ci1$percent[5], degree_ci2$percent[5]))

########
set.seed(2021)
degree_ci1 <- boot.ci(boot(betw_period[period == 1, betw], bootMean, R = 10000), type = "perc")
set.seed(2021)
degree_ci2 <- boot.ci(boot(betw_period[period == 2, betw], bootMean, R = 10000), type = "perc")

betw_error <- data.table(group = c("I", "II"),
                         value = c(degree_ci1$t0, degree_ci2$t0),
                         low = c(degree_ci1$percent[4], degree_ci2$percent[4]),
                         upper = c(degree_ci1$percent[5], degree_ci2$percent[5]))
betw_error
##########################

set.seed(2021)
degree_ci1 <- boot.ci(boot(diam_period[period == 1, dist], bootMean, R = 10000), type = "perc")
set.seed(2021)
degree_ci2 <- boot.ci(boot(diam_period[period == 2, dist], bootMean, R = 10000), type = "perc")

diam_error <- data.table(group = c("I", "II"),
                         value = c(degree_ci1$t0, degree_ci2$t0),
                         low = c(degree_ci1$percent[4], degree_ci2$percent[4]),
                         upper = c(degree_ci1$percent[5], degree_ci2$percent[5]))

######################
set.seed(2021)
degree_ci1 <- boot.ci(boot(size_period[period == 1, dist], bootMean, R = 10000), type = "perc")
set.seed(2021)
degree_ci2 <- boot.ci(boot(size_period[period == 2, dist], bootMean, R = 10000), type = "perc")

size_error <- data.table(group = c("I", "II"),
                         value = c(degree_ci1$t0, degree_ci2$t0),
                         low = c(degree_ci1$percent[4], degree_ci2$percent[4]),
                         upper = c(degree_ci1$percent[5], degree_ci2$percent[5]))


######## CI for totals
set.seed(2021)
boot.ci(boot(degree_period$out_degree, bootMean, R = 10000), type = "perc")
set.seed(2021)
boot.ci(boot(dist_period$dist, bootMean, R = 10000), type = "perc")
set.seed(2021)
boot.ci(boot(betw_period$betw, bootMean, R = 10000), type = "perc")
set.seed(2021)
boot.ci(boot(diam_period$dist, bootMean, R = 10000), type = "perc")
set.seed(2021)
boot.ci(boot(size_period$dist, bootMean, R = 10000), type = "perc")

######################

col <- c(pal_jama()(2)[1], viridis_pal(begin = 1, end = 0)(20)[15])
font_size <- 30
expr <- p_value_expr(test_results[measures == "avg_outdegree", p_value])
p_degree_error <- ggplot(degree_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    # scale_fill_jama() + 
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 0.97, size = 9, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_y_continuous(breaks = seq(0, 1.00, 0.25),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 1.00)) +
    labs(x = 'Period', y = 'Average out-degree') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = font_size, 
                                                        family = FONTS),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

expr <- p_value_expr(test_results[measures == "avg_distance", p_value])
p_dist_error <- ggplot(dist_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    # scale_fill_jama() + 
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 1.94, size = 9, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_y_continuous(breaks = seq(0, 2, 0.50),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 2)) +
    labs(x = 'Period', y = 'Average shortest path length') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = font_size, 
                                                        family = FONTS),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

expr <- p_value_expr(test_results[measures == "avg_betweenness", p_value])
p_betw_error <- ggplot(betw_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    # scale_fill_jama() + 
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 0.97, size = 9, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 1)) +
    labs(x = 'Period', y = 'Average betweenness') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = font_size, 
                                                        family = FONTS),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

expr <- p_value_expr(test_results[measures == "avg_diameter", p_value])
p_diam_error <- ggplot(diam_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    # scale_fill_jama() + 
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 1.552, size = 9, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_y_continuous(breaks = seq(0, 1.6, 0.40),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 1.6)) +
    labs(x = 'Period', y = 'Average diameter of clusters') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = font_size, 
                                                        family = FONTS),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

expr <- p_value_expr(test_results[measures == "avg_cluster_size", p_value])
p_size_error <- ggplot(size_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    # scale_fill_jama() + 
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 5.82, size = 9, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_y_continuous(breaks = seq(0, 6, 1.5),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 6)) +
    labs(x = 'Period', y = 'Average size of clusters') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = font_size, 
                                                        family = FONTS),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

# before and after network

out_degree <- sna::degree(net_before, cmode = 'freeman')
col <- c('#ffeda0', viridis_pal(begin = 1, end = 0)(20)[c(1, 5, 7, 10, 15, 20)])
names(col) <- as.character(c('0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '>=60'))
y <- plyr::mapvalues(network::get.vertex.attribute(net_before, 'age_interval'), from = 1:7,
                     to = c('0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '>=60'))

set.seed(2021)
g2 <- ggnet2(net_before, color = y, palette = col, color.legend = "Age group",
             size = 'degree', max_size = 8, arrow.gap = 0.005, arrow.size = 3, edge.color = "household",
             edge.size = 0.5,
             legend.size = 35) +
    guides(size = FALSE,
           color = guide_legend(size = 1, override.aes = list(size = 5))) +
    #ggtitle("Overall observed transmission network") +
    theme(text = element_text(family = FONTS),
          legend.position = 'left',
          legend.key.size = unit(1, "cm"),
          legend.key.height = unit(3, "line"))

out_degree <- sna::degree(net_after, cmode = 'freeman')
col <- c('#ffeda0', viridis_pal(begin = 1, end = 0)(20)[c(1, 5, 7, 10, 15, 20)])
names(col) <- as.character(c('0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '>=60'))
y <- plyr::mapvalues(network::get.vertex.attribute(net_after, 'age_interval'), from = 1:7,
                     to = c('0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '>=60'))

set.seed(2023)
g3 <- ggnet2(net_after, color = y, palette = col, color.legend = "Age group",
             size = 'degree', max_size = 6, arrow.gap = 0.006, arrow.size = 3, edge.color = "household",
             edge.size = 0.5,
             legend.size = 35) +
    guides(size = FALSE,
           color = guide_legend(size = 1, override.aes = list(size = 5))) +
    #ggtitle("Overall observed transmission network") +
    theme(text = element_text(family = FONTS),
          legend.position = 'left',
          legend.key.size = unit(1, "cm"),
          legend.key.height = unit(3, "line"))

Figure3 <- ggarrange(g2, g3, font.label = list(size = 30, face = 'plain', family=FONTS), widths = c(11, 11),
                     labels = c('(a)', '(b)'), common.legend = T, legend = "left")
Figure3
ggsave(filename = 'figure/Figure3.svg', width = 24, height = 11)
ggsave(filename = 'figure/Figure3.pdf', width = 24, height = 11, device=cairo_pdf)

Figure4 <- ggarrange(p_degree_error, p_dist_error, p_betw_error, 
                     p_diam_error, p_size_error, nrow = 1,
                     widths = c(8, 8, 8, 8, 8), labels = c('(a)', '(b)', '(c)', '(d)', '(e)'), 
                     vjust = 1.2, hjust = 0,
                     common.legend = T, legend = 'top',
                     font.label = list(size = 23, face = 'plain', family=FONTS))
Figure4

# Figure2 <- ggarrange(Figure2ab, Figure2cdef, nrow = 2,
#                      heights = c(11, 9))
# Figure2


ggsave(filename = 'figure/Figure4.svg', width = 18, height = 8)
ggsave(filename = 'figure/Figure4.pdf', width = 18, height = 8, device=cairo_pdf)


#######################Curve of infections of reconstruction#######################
all_contact_infect <- data.table(read_parquet('result1/all_contact_infect.parquet'))
day_information <- data.table(read_parquet('result1/day_information.parquet'))
id_information <- data.table(read_parquet('result1/id_information.parquet'))
# outbreak_date <- seq(as.Date('2020-01-08'), as.Date('2020-02-23'), 1)
# active_case_day <- unlist(lapply(outbreak_date, active_case))
# active_case_day <- data.table(outbreak_date, active_case = active_case_day)
# outbreak_date <- seq(as.Date('2020-01-08'), as.Date('2020-02-23'), 1)
# new_infection_day <- unlist(lapply(outbreak_date, new_infection))
# new_infection_day <- data.table(outbreak_date, new_infection_day)


###################################Dynamics of graphical measures by scenarios###################################

simulation_summary_graph <- list()
for (i in 1:7) {
    print(i)
    all_contact_infect_i <- all_contact_infect[scenario == i]
    day_information_i <- day_information[scenario == i]
    id_information_i <- id_information[scenario == i]
    Rt <- Rt_dt(i, 200, day_information_i, all_contact_infect_i, id_information_i, wuhan_date)
    simulation_summary_graph[[i]] <- Rt
}

simulation_summary_graph_total <- rbindlist(simulation_summary_graph)

simulation_summary_graph_python <- fread('result1/simulation_summary_total.csv')
simulation_summary_graph_python <- 
  simulation_summary_graph_python[!Measures %in% c("Number of singletons", "Percentage of singletons")]
simulation_summary_graph_total <- rbind(simulation_summary_graph_python, simulation_summary_graph_total)
simulation_summary_graph_total[, scenario := as.factor(scenario)]

write.csv(simulation_summary_graph_total, file = 'result1/simulation_summary_total_all.csv')
######################
simulation_summary_graph_total <- fread('result1/simulation_summary_total_all.csv')
simulation_summary_graph_total[, ymin := 0.00]
simulation_summary_graph_total[, ymax := 0.00]
simulation_summary_graph_total[Measures == 'Percentage of infection', ymin := 0]
simulation_summary_graph_total[Measures == 'Percentage of infection', ymax := 0.8]
simulation_summary_graph_total[Measures == 'Number of the new-onset', ymin := 0]
simulation_summary_graph_total[Measures == 'Number of the new-onset', ymax := 620]
simulation_summary_graph_total[Measures == 'Average outdegree', ymin := 0.5]
simulation_summary_graph_total[Measures == 'Average outdegree', ymax := 0.95]
simulation_summary_graph_total[Measures == 'Average diameter of clusters', ymin := 1]
simulation_summary_graph_total[Measures == 'Average diameter of clusters', ymax := 4]
simulation_summary_graph_total[Measures == 'Average size of clusters', ymin := 0]
simulation_summary_graph_total[Measures == 'Average size of clusters', ymax := 30]
simulation_summary_graph_total[Measures == 'Average shortest path length', ymin := 1]
simulation_summary_graph_total[Measures == 'Average shortest path length', ymax := 2.5]
simulation_summary_graph_total[Measures == 'Average betweenness', ymin := 0]
simulation_summary_graph_total[Measures == 'Average betweenness', ymax := 8]
simulation_summary_graph_total[Measures == 'Proportion within household', ymin := 0]
simulation_summary_graph_total[Measures == 'Proportion within household', ymax := 0.3]
simulation_summary_graph_total[Measures == 'Effective reproduction numbers', ymin := 0]
simulation_summary_graph_total[Measures == 'Effective reproduction numbers', ymax := 0.3]
######################
scenarios_label <- c('1: C(H) + L(S) + R(S) + HQ',
                     '2: C(H) + L(M) + R(M)',
                     '3: C(H) + L(M) + R(M) + HQ',
                     '4: C(L) + L(N) + R(M)',
                     '5: C(L) + L(M) + R(M)',
                     '6: C(L) + L(N) + R(M) + HQ',
                     '7: C(L) + L(M) + R(M) + HQ')
tmp <- simulation_summary_graph_total[Measures %in% c("Percentage of infection", 
                                                      "Number of the new-onset",
                                                      "Proportion within household",
                                                      "Effective reproduction numbers")]
tmp[, Measures := factor(Measures, levels = unique(Measures),
                         labels = c('A: Perentage of infection',
                                    'B: Number of the new-onset',
                                    'C: Proportion within household',
                                    'D: Effective reproduction numbers'))]
tmp[, scenario := factor(scenario, labels = scenarios_label)]
# col <- c("#B24745FF",
#          "#374E55FF", '#2ca25f',
#          "#e6550d", '#feb24c',
#          "#00A1D5FF", '#9ecae1')
# linetypes <- c('solid', 'solid', 'solid', 'dotted', 'dotted', 'dashed', 'dashed')

col <- c("black",
         "#2ca25f", '#2ca25f',
         "#e6550d", '#e6550d',
         "#00A1D5FF", '#00A1D5FF')
linetypes <- c('solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid')

ggplot(data = tmp, aes(color = scenario, linetype = scenario)) +
  geom_line(aes(x = day, y = median), size = 1.2) +
  geom_ribbon(aes(x=day, ymin=min, ymax=max, fill = scenario), alpha=.3) +
  annotate("rect", xmin = 16, xmax = 32, ymin = -Inf, ymax = Inf,
           alpha = .2) +
  facet_wrap_custom(~Measures, ncol = 2, scale = 'free_y',
                    scale_overrides = 
                      list(scale_override(1, scale_y_continuous(labels = scales::percent_format())), 
                           scale_override(3, scale_y_continuous(labels = scales::percent_format())), 
                           scale_override(2, scale_y_continuous(labels = seq(0, 400, 100),
                                                                breaks = seq(0, 400, 100), limits=c(0, 400))),
                           scale_override(4, scale_y_continuous(labels = seq(0, 5, 1),
                                                                breaks = seq(0, 5, 1), limits=c(0, 5))))) +
  # facet_grid_sc(rows = vars(Measures), cols=2, scales = list(y=scales_y))+
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  scale_linetype_manual(values = linetypes) +
  labs(color = 'Scenario', linetype = 'Scenario', x = "day") +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1.5),
                              nrow = 4, byrow = T, title = 'Scenarios'),
         fill = guide_legend(title = 'Scenarios'),
         linetype = guide_legend(title = 'Scenarios')) + 
  theme_bw() +
  theme(text = element_text(size = 26, family = FONTS),
        legend.title = element_text(size = 26, family = FONTS),
        axis.title=element_text(size = 26, family = FONTS),
        legend.position = 'top',
        axis.title.y = element_blank(),
        panel.border = element_rect(size = 1.1),
        legend.key.width = unit(3.5, "line"))
ggsave(filename = 'figure1/FigureS51.pdf', width = 12, height = 10, device=cairo_pdf)
######################
tmp <- simulation_summary_graph_total[Measures %in% c("Average outdegree", 
                                                      "Average shortest path length",
                                                      "Average betweenness",
                                                      "Average diameter of clusters",
                                                      "Average size of clusters")]
tmp[, Measures := factor(Measures, levels = unique(Measures),
                         labels = c('A: Average outdegree',
                                    'B: Average shortest path length',
                                    'C: Average betweenness',
                                    'D: Average diameter of clusters',
                                    'E: Average size of clusters'))]
tmp[, scenario := factor(scenario, labels = scenarios_label)]

ggplot(data = tmp, aes(color = scenario, linetype = scenario)) +
  geom_line(aes(x = day, y = median), size = 1.2) +
  geom_ribbon(aes(x=day, ymin=min, ymax=max, fill = scenario), alpha=.3) +
  annotate("rect", xmin = 16, xmax = 32, ymin = -Inf, ymax = Inf,
           alpha = .2) +
  facet_wrap_custom(~Measures, ncol = 2, scale = 'free_y',
                    scale_overrides = 
                      list(scale_override(1, scale_y_continuous(labels = seq(0.5, 1, 0.1),
                                                                breaks = seq(0.5, 1, 0.1),
                                                                limits = c(0.5, 1))),
                           scale_override(2, scale_y_continuous(labels = seq(1, 5, 1),
                                                                breaks = seq(1, 5, 1),
                                                                limits = c(1, 5))),
                           scale_override(3, scale_y_continuous(labels = seq(0, 21, 3),
                                                                breaks = seq(0, 21, 3),
                                                                limits = c(0, 21))),
                           scale_override(4, scale_y_continuous(labels = seq(1, 6, 1),
                                                                breaks = seq(1, 6, 1),
                                                                limits = c(1, 6.5))),
                           scale_override(5, scale_y_continuous(labels = seq(0, 40, 10),
                                                                breaks = seq(0, 40, 10),
                                                                limits = c(0, 40))))) +
  # facet_grid_sc(rows = vars(Measures), cols=2, scales = list(y=scales_y))+
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  scale_linetype_manual(values = linetypes) +
  labs(color = 'Scenario', linetype = 'Scenario', x = "day") +
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1.5),
                              nrow = 4, byrow = T, title = 'Scenarios'),
         fill = guide_legend(title = 'Scenarios'),
         linetype = guide_legend(title = 'Scenarios')) + 
  theme_bw() +
  theme(text = element_text(size = 26, family = FONTS),
        legend.title = element_text(size = 26, family = FONTS),
        axis.title=element_text(size = 26, family = FONTS),
        legend.position = 'top',
        axis.title.y = element_blank(),
        panel.border = element_rect(size = 1.1),
        legend.key.width = unit(3.5, "line"))
ggsave(filename = 'figure1/FigureS52.pdf', width = 12, height = 14, device=cairo_pdf)

# 
# ggplot(data = tmp, aes(color = scenario, linetype = scenario)) +
#     geom_line(aes(x = day, y = median), size = 1.2) +
#     # geom_ribbon(aes(x=day, ymin=min, ymax=max, fill = scenario), alpha=.3) +
#     annotate("rect", xmin = 16, xmax = 32, ymin = -Inf, ymax = Inf,
#              alpha = .2) +
#     facet_wrap_custom(~Measures, ncol = 2, scale = 'free_y',
#                       scale_overrides = list(scale_override(1,
#                                                             scale_y_continuous(labels = scales::percent_format())),
#                                              scale_override(7,
#                                                             scale_y_continuous(labels = scales::percent_format())),
#                                              scale_override(2,
#                                                             scale_y_continuous(labels = seq(0, 600, 200),
#                                                                                breaks = seq(0, 600, 200),
#                                                             limits=c(0, 600))))) +
#     facet_grid_sc(rows = vars(Measures), cols=2, scales = list(y=scales_y))+
#     geom_blank(aes(y = ymin)) +
#     geom_blank(aes(y = ymax)) +
#     scale_color_manual(values = col) +
#     scale_linetype_manual(values = linetypes) +
#     labs(color = 'Scenario', linetype = 'Scenario', x = "天数") +
#     guides(color = guide_legend(override.aes = list(size = 1, alpha = 1.5),
#                                 nrow = 4, byrow = T, title = '仿真设定'),
#            # fill = guide_legend(title = '仿真设定'),
#            linetype = guide_legend(title = '仿真设定')) + 
#     theme_bw() +
#     theme(text = element_text(size = 24, family = FONTS),
#           legend.title = element_text(size = 20, family = FONTS),
#           axis.title=element_text(size=20, family = FONTS),
#           legend.position = 'top',
#           axis.title.y = element_blank(),
#           panel.border = element_rect(size = 1.1),
#           legend.key.width = unit(3.5, "line"))
#     
# ggsave(filename = 'figure/Figure5.svg', width = 12, height = 14)
# ggsave(filename = 'figure/Figure5.pdf', width = 12, height = 16, device=cairo_pdf)


###################################Social description###################################
social_summary <- data.table(read_parquet('/Users/howanchung/Downloads/net_attr_day.parquet'))
# social_summary <- fread('/Users/howanchung/Downloads/social_summary_total.csv')
setnames(social_summary, colnames(social_summary), 
         c("Average outdegree", "Average shortest path length", 
           "Average diameter of clusters",  
           "Proportion within household",
           "seed", "scenario", "day"))

social_summary <- melt(social_summary, id.vars = c("day", "scenario", "seed"), 
     measure.vars = c("Average outdegree", "Average shortest path length", 
                      "Average diameter of clusters", "Proportion within household"),
     variable.name = "Measures", value.name = "x")

social_summary <- social_summary %>%
  group_by(day, scenario, Measures) %>%
  summarise(min=quantile(x, 0.025),
            median = median(x),
            max=quantile(x, 0.975)) %>%
  as.data.table()

social_summary <- social_summary[Measures %in% 
                                   c('Average outdegree', 'Average diameter of clusters',
                                     'Average shortest path length',
                                     'Proportion within household'),]
social_summary <- social_summary[day <= 47,]
social_summary[, scenario := factor(scenario)]
# social_summary[, Measures := factor(Measures,
#                                     levels = unique(social_summary$Measures),
#                                     labels = c('A: Average outdegree',
#                                                'B: Average diameter of clusters',
#                                                'C: Average size of clusters',
#                                                'D: Average shortest path length',
#                                                'E: Proportion within household'))]


simulate_dynamics <- ggplot(data = social_summary, aes(color = scenario, linetype = scenario)) +
  geom_line(aes(x = day, y = median), size = 0.9) +
  geom_ribbon(aes(x=day, ymin=min, ymax=max, fill = scenario), alpha=.3) + 
  facet_wrap_custom(~Measures, ncol = 2, scale = 'free_y',
                    scale_overrides = list(
                      scale_override(4, scale_y_continuous(labels = scales::percent_format())))) +
  scale_fill_manual(values = c("#d95f0e", "#2c7fb8")) +
  theme_bw() +
  labs(color = 'Scenario', linetype = 'Scenario', fill = 'Scenario') +
  theme(text = element_text(size = 25),
        legend.position = 'top',
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(size = 1.1)) +
  guides(color = guide_legend(override.aes = list(size = 0.7),
                              nrow = 2, byrow = T))

simulate_dynamics



