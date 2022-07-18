summary_graph <- function(T_network, feature = 'diameter') {
    singletons <- which(sna::degree(T_network) == 0)
    T_graph <- intergraph::asIgraph(T_network)
    # T_graph <- igraph::delete_vertices(T_graph, singletons)
    if (feature == 'diameter') {
        return(sapply(igraph::decompose(T_graph), igraph::diameter))
    }else if (feature == 'cluster size') {
        return(sapply(igraph::decompose(T_graph), igraph::vcount))
    }else if (feature == 'distance') {
        distances <- igraph::distance_table(T_graph, T)$res
        return(rep(1:length(distances), distances))
    }else {
        return(FALSE)
    }
}


in_degree <- sna::degree(net, cmode = 'indegree')
out_degree <- sna::degree(net, cmode = 'outdegree')
nonsingletons <- 1:1350
degree_period <- data.table(id = 1:1350, out_degree = out_degree)
degree_period[, period := 2]
degree_period[id %in% (1:1350)[COVID19$origin_time <= '2020-01-23'], period := 1]
degree_period <- degree_period[id %in% nonsingletons]

degree_period %>%
    group_by(period) %>%
    summarise(mean_outdegree = mean(out_degree))


########
degree_dist <- degree_period %>%
    group_by(period) %>%
    summarise('0' = sum(out_degree == 0),
              '1' = sum(out_degree == 1),
              '2' = sum(out_degree == 2),
              '>=3' = sum(out_degree >= 3)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "degree") %>%
    as.data.table()
degree_dist[, total := sum(N), by = period][, prop := N / total][, total := NULL]
degree_dist[, period := factor(period, levels = c(1, 2), labels = c("I", "II"))]

set.seed(2021)
degree_ci1 <- boot.ci(boot(degree_period[period == 1, out_degree], bootMean, R = 10000), type = "norm")
set.seed(2021)
degree_ci2 <- boot.ci(boot(degree_period[period == 2, out_degree], bootMean, R = 10000), type = "norm")

degree_error <- data.table(group = c("I", "II"),
                           value = c(degree_ci1$t0, degree_ci2$t0),
                           low = c(degree_ci1$normal[2], degree_ci2$normal[2]),
                           upper = c(degree_ci1$normal[3], degree_ci2$normal[3]))

########

dist_period <- rbind(data.table(period = 1, dist = summary_graph(net_before, feature = 'distance')),
                     data.table(period = 2, dist = summary_graph(net_after, feature = 'distance')))

dist_dist <- dist_period %>%
    group_by(period) %>%
    summarise('1' = sum(dist == 1),
              '2' = sum(dist == 2),
              '3' = sum(dist == 3),
              '4' = sum(dist == 4)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "dist") %>%
    as.data.table()
dist_dist[, total := sum(N), by = period][, prop := N / total][, total := NULL]
dist_dist[, period := factor(period, levels = c(1, 2), labels = c("I", "II"))]

set.seed(2021)
degree_ci1 <- boot.ci(boot(dist_period[period == 1, dist], bootMean, R = 10000), type = "norm")
set.seed(2021)
degree_ci2 <- boot.ci(boot(dist_period[period == 2, dist], bootMean, R = 10000), type = "norm")

dist_error <- data.table(group = c("I", "II"),
                         value = c(degree_ci1$t0, degree_ci2$t0),
                         low = c(degree_ci1$normal[2], degree_ci2$normal[2]),
                         upper = c(degree_ci1$normal[3], degree_ci2$normal[3]))

##########################

diam_period <- rbind(data.table(period = 1, dist = summary_graph(net_before, feature = 'diameter')),
                     data.table(period = 2, dist = summary_graph(net_after, feature = 'diameter')))

diam_dist <- diam_period %>%
    group_by(period) %>%
    summarise('1' = sum(dist == 1),
              '2' = sum(dist == 2),
              '3' = sum(dist == 3),
              '4' = sum(dist == 4)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "diam") %>%
    as.data.table()
diam_dist[, total := sum(N), by = period][, prop := N / total][, total := NULL]
diam_dist[, period := factor(period, levels = c(1, 2), labels = c("I", "II"))]

set.seed(2021)
degree_ci1 <- boot.ci(boot(diam_period[period == 1, dist], bootMean, R = 10000), type = "norm")
set.seed(2021)
degree_ci2 <- boot.ci(boot(diam_period[period == 2, dist], bootMean, R = 10000), type = "norm")

diam_error <- data.table(group = c("I", "II"),
                         value = c(degree_ci1$t0, degree_ci2$t0),
                         low = c(degree_ci1$normal[2], degree_ci2$normal[2]),
                         upper = c(degree_ci1$normal[3], degree_ci2$normal[3]))

######################

size_period <- rbind(data.table(period = 1, dist = summary_graph(net_before, feature = 'cluster size')),
                     data.table(period = 2, dist = summary_graph(net_after, feature = 'cluster size')))

size_dist <- size_period %>%
    group_by(period) %>%
    summarise('2' = sum(dist == 2),
              '3' = sum(dist == 3),
              '4' = sum(dist == 4),
              '>=5' = sum(dist >= 5)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "size") %>%
    as.data.table()
size_dist[, total := sum(N), by = period][, prop := N / total][, total := NULL]
size_dist[, period := factor(period, levels = c(1, 2), labels = c("I", "II"))]

set.seed(2021)
degree_ci1 <- boot.ci(boot(size_period[period == 1, dist], bootMean, R = 10000), type = "norm")
set.seed(2021)
degree_ci2 <- boot.ci(boot(size_period[period == 2, dist], bootMean, R = 10000), type = "norm")

size_error <- data.table(group = c("I", "II"),
                         value = c(degree_ci1$t0, degree_ci2$t0),
                         low = c(degree_ci1$normal[2], degree_ci2$normal[2]),
                         upper = c(degree_ci1$normal[3], degree_ci2$normal[3]))


######## CI for totals
set.seed(2021)
boot.ci(boot(degree_period$out_degree, bootMean, R = 10000), type = "norm")
set.seed(2021)
boot.ci(boot(dist_period$dist, bootMean, R = 10000), type = "norm")
set.seed(2021)
boot.ci(boot(diam_period$dist, bootMean, R = 10000), type = "norm")
set.seed(2021)
boot.ci(boot(size_period$dist, bootMean, R = 10000), type = "norm")

######## t-test
t.test(out_degree ~ period, data = degree_period, alternative = 'greater')
t.test(x = summary_graph(net_before, feature = 'diameter'),
       y = summary_graph(net_after, feature = 'diameter'),
       alternative = 'greater')
t.test(x = summary_graph(net_before, feature = 'distance'),
       y = summary_graph(net_after, feature = 'distance'),
       alternative = 'greater')
t.test(x = summary_graph(net_before, feature = 'cluster size'),
       y = summary_graph(net_after, feature = 'cluster size'),
       alternative = 'greater')

######################

col <- c(pal_jama()(2)[1], viridis_pal(begin = 1, end = 0)(20)[15])

p_degree_error <- ggplot(degree_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 0.57, size = 9, label = "p-value < 0.001",
             fontface = 'italic') +
    scale_y_continuous(breaks = seq(0, 0.6, 0.15),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 0.6)) +
    labs(x = 'Period', y = 'Average out-degree') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 25),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

p_dist_error <- ggplot(dist_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 1.90, size = 9, label = "p-value < 0.001",
             fontface = 'italic') +
    scale_y_continuous(breaks = seq(0, 2, 0.50),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 2)) +
    labs(x = 'Period', y = 'Average shortest path length') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 25),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

p_diam_error <- ggplot(diam_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 0.57, size = 9, label = "p-value < 0.001",
             fontface = 'italic') +
    scale_y_continuous(breaks = seq(0, 0.6, 0.15),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 0.6)) +
    labs(x = 'Period', y = 'Average diameter of clusters') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 25),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

p_size_error <- ggplot(size_error) +
    geom_bar(aes(x = fct_inorder(group), y = value, fill = group), stat = 'identity',
             width = .75) +
    scale_fill_manual(values = col) +
    geom_errorbar(aes(x = group, ymin = low, ymax = upper), width = .2,
                  position = position_dodge(.9)) +
    annotate('text', x = 1.5, y = 2.28, size = 9, label = "p-value < 0.001",
             fontface = 'italic') +
    scale_y_continuous(breaks = seq(0, 2.4, 0.6),
                       labels = scales::number_format(accuracy = 0.01),
                       limits = c(0, 2.4)) +
    labs(x = 'Period', y = 'Average sizes of clusters') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 25),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

Figure4 <- ggarrange(p_degree_error, p_dist_error, p_diam_error, p_size_error, nrow = 1,
                     widths = c(8, 8, 8, 8), labels = c('(a)', '(b)', '(c)', '(d)'), vjust = 1.2, hjust = 0,
                     common.legend = T, legend = 'top',
                     font.label = list(size = 23, face = 'plain'))
Figure4

ggsave(filename = 'result/FigureS10.svg', width = 15, height = 8)
ggsave(filename = 'result/FigureS10.pdf', width = 18, height = 8)