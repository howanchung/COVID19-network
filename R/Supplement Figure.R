source("R/data preparation.R")
FONTS <- "STSong"
########################

age_suscep <- data.table(age_group = fct_inorder(c('0-4', '5-9', '10-14', '15-19',
                                                   '20-24', '25-29', '30-34', '35-39',
                                                   '40-44', '45-49', '50-54', '55-59',
                                                   '60-64', '>=65')),
                         suscep = c(rep(c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82), each = 2), 0.88, 0.81))

col <- pal_jama()(2)[1]
ggplot(data = age_suscep, aes(x = age_group, y = suscep)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    geom_text(aes(x = age_group, y = suscep, label = suscep), 
              color = 'black', size = 8, vjust = -0.5,
              family = FONTS) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) + 
    labs(x = 'Age group', y = 'Susceptibility') +
    theme(text = element_text(size = 30, family = FONTS),
          axis.text=element_text(size = 26),
          panel.border = element_rect(size = 1.1),
          legend.position = "none")

ggsave(filename = 'figure/FigureS1.svg', width = 16, height = 8)
ggsave(filename = 'figure/FigureS1.pdf', width = 14, height = 7, device = cairo_pdf)

########################
contact_curve <- function(day, begin, baseline, eta, lambda, m) {
    baseline * ((1 - eta) / (1 + exp(lambda * (day - begin - m / 2))) + eta)
}

mean_contact_curve <- data.table(day = 1:20,
                                 contact_rate = sapply(1:20, contact_curve, begin = 3, baseline = 19.176,
                                                       eta = 2.215 / 19.176, lambda = 2 * log(99) / 13, m = 13))

ggplot(data = mean_contact_curve) +
    geom_line(aes(x = day, y = contact_rate)) +
    scale_y_continuous(breaks = c(mean_contact_curve[1, contact_rate], mean_contact_curve[16, contact_rate]),
                       labels = c(TeX(r'($\textit{c_{ij}^{(0)}}$)'), TeX(r'($\textit{\eta_{ij}c_{ij}^{(0)}}$)'))) +
    scale_x_continuous(breaks = c(3, 16), labels = c(TeX(r'($\textit{{t_0}}$)'), TeX(r'($\textit{t_0+m}$)'))) +
    theme_bw() +
    labs(x = 'Day', y = 'Number of social contacts') +
    theme(text = element_text(size = 30),
          axis.text=element_text(size = 28),
          axis.title.x = element_text(size=30, family = FONTS),
          axis.title.y = element_text(size=30, family = FONTS))

ggsave(filename = 'figure/FigureS2.svg', width = 12, height = 8)
ggsave(filename = 'figure/FigureS2.pdf', width = 14, height = 8, device = cairo_pdf)


#########################
contact_pattern <- fread("result/contact_pattern.csv")
setnames(contact_pattern, c('0-4', '5-9', '10-14', '15-19',
                            '20-24', '25-29', '30-34', '35-39',
                            '40-44', '45-49', '50-54', '55-59',
                            '60-64', '>=65'))
contact_pattern[, day := seq(0, nrow(contact_pattern) - 1) + as.Date('2020-01-08')]
contact_pattern <- contact_pattern[day <= (as.Date('2020-01-08') + 64 + 5)]
contact_pattern <- melt(data = contact_pattern, id.vars = "day",
                        value.name = 'contact_rate', variable.name = "age_group")

ggplot(data = contact_pattern, aes(x = day)) +
    geom_line(aes(x = day, y = contact_rate, color = age_group)) +
    annotate("rect", xmin = as.Date('2020-01-08') + 16,
             xmax = as.Date('2020-01-08') + 33, ymin = 0, ymax = Inf, alpha = .2) +
    annotate('text', x = as.Date('2020-01-08') + 25,
             y = 25, size = 8, label = "Period with the highest-level alert", family = FONTS) +
    scale_x_date(breaks = c(0, 16, 33, 63) + as.Date('2020-01-08'), date_labels = "%b %d") +
    labs(x = 'Day', y = 'Average number of social contacts', color = 'Age group') +
    theme_bw() +
    theme(text = element_text(size = 30, family = FONTS),
          legend.position = 'top',
          panel.border = element_rect(size = 1.1)) +
    scale_color_simpsons() +
    guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = 'figure/FigureS3.svg', width = 16, height = 8)
ggsave(filename = 'figure/FigureS3.pdf', width = 14, height = 7, device = cairo_pdf)

########################total number of cases by age groups
age_cases <- COVID19[, c("age_year")]
age_cases[, age_group := cut(as.numeric(age_year), breaks = c(seq(0, 60, 10), 100), right = FALSE, labels = FALSE)]
age_cases[, age_group := factor(age_group, levels = 1:7, labels = AGE_LABEL)]

age_cases <- age_cases %>%
  group_by(age_group) %>%
  summarise(count = n()) %>%
  as.data.table()

figure_age_cases <- ggplot(age_cases, aes(x = fct_inorder(age_group), y = count)) +
  geom_bar(stat = 'identity', width = .75) +
  geom_text(aes(label = count, y = count),
            vjust = -0.5, size = 7, family = FONTS) +
  scale_y_continuous(breaks = seq(0, 350, 70),
                     limits = c(0, 350)) +
  labs(x = 'Age group', y = 'Number of cases') +
  theme_bw() +
  theme(legend.position = 'left', text = element_text(size = 30, family = FONTS),
        panel.border = element_rect(size = 1.1))
figure_age_cases
ggarrange(figure_age_cases, font.label = list(size = 24, face = 'plain', family = FONTS),
          labels = '(a)', common.legend = TRUE, vjust=1.5)
ggsave(filename = 'figure/age_cases.pdf', width = 14, height = 7, device = cairo_pdf)

#########################

age_infect <- data.table(aggregate(COVID19_CJ$infect, by = list(age_group1 = COVID19_CJ$age_group1,
                                                                age_group2 = COVID19_CJ$age_group2),
                                   FUN = sum))
age_infect[, age_group1 := factor(age_group1, levels = 1:7, labels = AGE_LABEL)]
age_infect[, age_group2 := factor(age_group2, levels = 1:7, labels = AGE_LABEL)]

FigureS5a <- ggplot(data = age_infect, aes(x = fct_inorder(age_group2), y = fct_inorder(age_group1), fill = x)) +
    geom_tile() +
    scale_fill_gradient(low = 'white', high = 'red') +
    theme_bw() +
    geom_text(aes(age_group2, age_group1, label = x), color = 'black', size = 6, family = FONTS) +
    labs(x = 'Age group of secondary cases', y = 'Age group of source cases', fill = NULL) +
    theme(text = element_text(size = 30, family = FONTS),
          panel.border = element_blank(),
          legend.position = "none")
FigureS5a

ggsave(filename = 'figure/age_age.pdf', width = 12, height = 9, device = cairo_pdf)


############
in_degree <- sna::degree(net, cmode = 'indegree')
out_degree <- sna::degree(net, cmode = 'outdegree')
nonsingletons <- which(in_degree != 0 | out_degree != 0)
outdegree_age <- data.table(id = nonsingletons)
outdegree_age[, age := plyr::mapvalues(id,
                                       from = COVID19_CJ$id1,
                                       to = COVID19_CJ$age_group1, warn_missing = F)]
outdegree_age[, outdegree := out_degree[nonsingletons]]

# CI for out-degree for cases aged from 40 to 49
mean(outdegree_age[age == 5, outdegree])
set.seed(2021)
boot.ci(boot(outdegree_age[age == 5, outdegree], bootMean, R = 10000), type = "perc")

outdegree_age <- outdegree_age %>%
    group_by(age) %>%
    summarise(mean_outdegree = mean(outdegree),
              std = sd(outdegree) / sqrt(n()),
              n = n()) %>%
    as.data.table()


outdegree_age[, age := factor(age, labels = AGE_LABEL)]

FigureS5b <- ggplot(outdegree_age, aes(x = fct_inorder(age), y = mean_outdegree)) +
    geom_bar(stat = 'identity', width = .75) +
    geom_errorbar(aes(x = fct_inorder(age), 
                    ymin = mean_outdegree - 1.96 * std, 
                    ymax = mean_outdegree + 1.96 * std), width = .2,
                position = position_dodge(.9)) +
    geom_text(aes(label = round(mean_outdegree, 2), y = mean_outdegree),
              vjust = -0.5, size = 7, family = FONTS) +
    scale_y_continuous(breaks = seq(0, 1.4, 0.2),
                       limits = c(0, 1.4)) +
    labs(x = 'Age group', y = 'Average out-degree') +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1))

############
between <- sna::betweenness(net)
betw_age <- data.table(id = nonsingletons)
betw_age[, age := plyr::mapvalues(id,
                                  from = COVID19_CJ$id1,
                                  to = COVID19_CJ$age_group1, warn_missing = F)]
betw_age[, betw_age := between[nonsingletons]]

betw_age <- betw_age %>%
    group_by(age) %>%
    summarise(avg_betw = mean(betw_age),
              std = sd(betw_age) / sqrt(n()),
              n = n()) %>%
    as.data.table()


betw_age[, age := factor(age, labels = AGE_LABEL)]

FigureS5c <- ggplot(betw_age, aes(x = fct_inorder(age), y = avg_betw)) +
    geom_bar(stat = 'identity', width = .75) +
    geom_errorbar(aes(x = fct_inorder(age), 
                      ymin = pmax(avg_betw - 1.96 * std, 0), 
                      ymax = avg_betw + 1.96 * std), width = .2,
                  position = position_dodge(.9)) +
    geom_text(aes(label = round(avg_betw, 2), y = avg_betw),
              vjust = -0.5, size = 6, family = FONTS) +
    scale_y_continuous(breaks = seq(0, 1.4, 0.2),
                       limits = c(0, 1.4)) +
    labs(x = 'Age group', y = 'Average betweeness centrality') +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1))

ggarrange(FigureS5b, FigureS5c, nrow=1,
          font.label = list(size = 24, face = 'plain', family = FONTS),
          labels = c('(b)', '(c)'), legend = "none", vjust=1.2)
ggsave(filename = 'figure/FigureS5ab.svg', width = 16, height = 8)
ggsave(filename = 'figure/FigureS5ab.pdf', width = 14, height = 7, device = cairo_pdf)

##
in_degree <- sna::degree(net, cmode = 'indegree')
out_degree <- sna::degree(net, cmode = 'outdegree')
nonsingletons <- which(in_degree != 0 | out_degree != 0)
nonsingletons_age <- data.table(id = nonsingletons)
nonsingletons_age[, age := plyr::mapvalues(id,
                                           from = COVID19_CJ$id1,
                                           to = COVID19_CJ$age_group1, warn_missing = F)]
indexed <- which(in_degree == 0 & out_degree > 0)
last_generation <- which(in_degree != 0 & out_degree == 0)

nonsingletons_age[, is_first := FALSE]
nonsingletons_age[id %in% indexed, is_first := TRUE]
nonsingletons_age[, is_last := FALSE]
nonsingletons_age[id %in% last_generation, is_last := TRUE]

# CI for proportion to being the indexed cases for cases aged from 40 to 49
mean(nonsingletons_age[age == 5, is_first])
set.seed(2021)
boot.ci(boot(nonsingletons_age[age == 5, is_first], bootMean, R = 10000), type = "perc")

age_index_last <- nonsingletons_age %>%
    group_by(age) %>%
    summarise(n_first = sum(is_first),
              prop_first = mean(is_first),
              se_first = sd(is_first) / sqrt(n()),
              n_last = sum(is_last),
              prop_last = mean(is_last),
              se_last = sd(is_last) / sqrt(n()),
              .groups = 'drop') %>%
    as.data.table()
age_index_last

nonsingletons_age %>%
    subset(!age %in% c(5, 6)) %>%
    summarise(n_last = sum(is_last),
              prop_last = mean(is_last),
              se_last = sd(is_last) / sqrt(n()))

# CI for proportion to being the indexed cases for cases aged from 40 to 49
mean(nonsingletons_age[!age %in% c(5, 6), is_last])
set.seed(2021)
boot.ci(boot(nonsingletons_age[!age %in% c(5, 6), is_last], bootMean, R = 10000), type = "perc")

nonsingletons_age <- as.data.table(table(nonsingletons_age$age))

generation_age <- age_index_last[, c("age", "n_last", "prop_last", "se_last")]
setnames(generation_age, c("age_group", "n", "prop", "se"))
index_age <- age_index_last[, c("age", "n_first", "prop_first", "se_first")]
setnames(index_age, c("age_group", "n", "prop", "se"))
age_index_last <- rbind(generation_age, index_age)
age_index_last[, age_group := factor(age_group, levels = 1:7,
                                     labels = AGE_LABEL)]
age_index_last[, type := c(rep('Terminal cases', 7), rep('Indexed cases', 7))]
age_index_last[, low := pmax(prop - 1.96 * se, 0)]
age_index_last[, upp := pmin(prop + 1.96 * se, 1)]

col <- pal_npg("nrc", alpha=0.8)(9)[c(4, 9)]
FigureS5d <- ggplot(age_index_last, aes(x = fct_inorder(age_group), y = prop, fill = type)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 5, family = FONTS) +
    geom_errorbar(aes(x = fct_inorder(age_group), ymin = low, ymax = upp), width = .2,
                  position = position_dodge(.9)) +
    # annotate('text', x = 1.75, y = 0.95, size = 10, family = 'Calibri', label = "p-value = 0.038",
    #          fontface = 'italic') +
    scale_fill_manual(values = col) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(x = 'Age group', y = 'Proportion') +
    guides(fill = guide_legend(title = 'Location in network')) +
    theme_bw() +
    theme(text = element_text(size = 30, family = FONTS),
          legend.position = 'top')

ggarrange(FigureS5d, font.label = list(size = 24, face = 'plain', family = FONTS),
          labels = '(d)', common.legend = TRUE, vjust=-0.5)
# FigureS5c
ggsave(filename = 'figure/FigureS5c.svg', width = 15, height = 8)
ggsave(filename = 'figure/FigureS5c.pdf', width = 14, height = 7, device = cairo_pdf)

##########
in_degree <- sna::degree(net, cmode = 'indegree')
out_degree <- sna::degree(net, cmode = 'outdegree')
trans_net3 <- COVID19_CJ[infect == 1, c('id1', 'id2', 'generation', 'household')]
trans_net3[is.na(household), household := 0]
trans_net3[, max_generation := plyr::mapvalues(trans_net3$id2, from = distance_dt$id2,
                                               to = distance_dt$max_generation, warn_missing = T)]
trans_net3[, max_generation := as.integer(max_generation)]
# trans_net3 <- trans_net3[generation >0, ]

indexed <- which(in_degree == 0 & out_degree > 0)
last_generation <- which(in_degree != 0 & out_degree == 0)
trans_net3[, is_first := (id1 %in% indexed)]
trans_net3[, is_last := (id2 %in% last_generation)]
trans_net3[, age_group1 := plyr::mapvalues(id1, from = COVID19_CJ$id1, to = COVID19_CJ$age_group1,
                                           warn_missing = F)]
trans_net3[, age_group2 := plyr::mapvalues(id2, from = COVID19_CJ$id2, to = COVID19_CJ$age_group2,
                                           warn_missing = F)]

trans_net3[, age_group1 := factor(age_group1, levels = 1:7,
                                  labels = AGE_LABEL)]
trans_net3[, age_group2 := factor(age_group2, levels = 1:7,
                                  labels = AGE_LABEL)]


# CI for proportion to household transmission
boot.ci(boot(trans_net3[, household], bootMean, R = 10000), type = "perc")
mean(trans_net3[, household])

# CI for proportion to being the terminal transmission for household transmission
boot.ci(boot(trans_net3[is_last == F, household], bootMean, R = 10000), type = "perc")
mean(trans_net3[is_last == T, household])

household_dist <- trans_net3 %>%
    group_by(is_last) %>%
    dplyr::summarise(ratio = mean(household),
                     se_ratio = sd(household) / sqrt(n()),
                     n = n()) %>%
    as.data.table()

household_dist[, is_last := c("FALSE", "TRUE")]

FigureS6a <- ggplot(household_dist, aes(x = is_last, y = ratio)) +
    geom_bar(stat = 'identity', width = .75) +
    geom_text(aes(label = scales::percent(ratio, accuracy = 0.1), y = ratio),
              vjust = -0.5, size = 6, family = FONTS) +
    geom_errorbar(aes(x = fct_inorder(is_last), 
                      ymin = ratio - 1.96 * se_ratio, 
                      ymax = ratio + 1.96 * se_ratio), width = .2,
                  position = position_dodge(.9)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(x = 'Terminal transmission', y = 'Proportion within household') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1))
FigureS6a


###################
household_dist2 <- trans_net3 %>%
    subset(max_generation > 1 & generation > 0) %>%
    group_by(generation) %>%
    dplyr::summarise(ratio = mean(household),
                     se_ratio = sd(household) / sqrt(n()),
                     n = n()) %>%
    as.data.table()

FigureS6b <- ggplot(household_dist2, aes(x = generation, y = ratio)) +
    geom_bar(stat = 'identity', width = .75) +
    geom_point(data = household_dist2, aes(x = generation, y = n / 200, group = 1), 
               color = pal_jama()(1)[1],
               size = 3, stat = 'identity') +
    geom_errorbar(aes(x = generation, 
                      ymin = pmax(ratio - 1.96 * se_ratio, 0), 
                      ymax = pmin(ratio + 1.96 * se_ratio, 1)), width = .2,
                  position = position_dodge(.9)) +
    geom_line(data = household_dist2, aes(x = generation, y = n / 200, group = 1),
              size = 1.5, stat = "identity") +
    geom_text(data = household_dist2,
              aes(x = generation, y = n / 200, label = n), vjust = -1, color = 'black', 
              size = 6, family = FONTS) +
    geom_text(aes(label = scales::percent(ratio, accuracy = 0.1), y = ratio),
              vjust = -0.5, size = 6, family = FONTS) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1),
                       sec.axis = sec_axis(~. * 200, name = 'Number of transmission')) +
    labs(x = 'Generation of transmission', y = 'Proportion within household') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'left', text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1))

ggarrange(FigureS6a, FigureS6b, font.label = list(size = 25, face = 'plain', family = FONTS),
          labels = c('(a)', '(b)'), legend = "none", widths = c(6, 10), vjust=1.2)
ggsave(filename = 'figure/FigureS6.svg', width = 16, height = 10)
ggsave(filename = 'figure/FigureS6.pdf', width = 14, height = 7, device = cairo_pdf)

############dynamic
########################################

trans_before_id <- COVID19_CJ_before %>%
    subset(infect != 0, select = id2)
trans_after_id <- COVID19_CJ_after %>%
    subset(infect != 0, select = id2)

trans_before <- trans_net3[id2 %in% plyr::mapvalues(trans_before_id$id2, from = COVID19$id2, to = 1:1349, warn_missing = F)]
trans_before[, period := 'I']
trans_after <- trans_net3[id2 %in% plyr::mapvalues(trans_after_id$id2, from = COVID19$id2, to = 1:1349, warn_missing = F)]
trans_after[, period := 'II']

trans_period <- rbind(trans_before, trans_after)

generation_age_period <- trans_period %>%
    group_by(is_last, age_group2, period) %>%
    summarise(n = n(), .groups = 'drop') %>%
    as.data.table() %>%
    subset(is_last == TRUE)

generation_age_period[, is_last := NULL]
generation_age_period[, type := rep('Terminal cases', 14)]
setnames(generation_age_period, c('age_group', 'Period', 'n', 'type'))


in_degree <- sna::degree(net_before, cmode = 'indegree')
out_degree <- sna::degree(net_before, cmode = 'outdegree')
indexed1 <- which(in_degree == 0 & out_degree > 0)
index_age1 <- data.table(id = network::get.vertex.attribute(net_before, 'vertex.names')[indexed1])
index_age1[, age_group := plyr::mapvalues(id, from = COVID19_CJ_before$id1,
                                          to = COVID19_CJ_before$age_group1, warn_missing = F)]
index_age1 <- as.data.table(table(index_age1$age_group))
index_age1 <- rbind(index_age1, data.table(V1 = 1, N = 0))
index_age1[, type := rep('Indexed cases', 7)]
setkey(index_age1, V1)
setnames(index_age1, c('age_group', 'n', 'type'))
index_age1[, age_group := factor(age_group, levels = 1:7,
                                 labels = c('0-9',
                                            '10-19',
                                            '20-29',
                                            '30-39',
                                            '40-49',
                                            '50-59',
                                            '>=60'))]

in_degree <- sna::degree(net_after, cmode = 'indegree')
out_degree <- sna::degree(net_after, cmode = 'outdegree')
indexed2 <- which(in_degree == 0 & out_degree > 0)
index_age2 <- data.table(id = network::get.vertex.attribute(net_after, 'vertex.names')[indexed2])
index_age2[, age_group := plyr::mapvalues(id, from = COVID19_CJ_after$id1,
                                          to = COVID19_CJ_after$age_group1, warn_missing = F)]
index_age2 <- as.data.table(table(index_age2$age_group))
index_age2 <- rbind(index_age2, data.table(V1 = 1, N = 0))
index_age2[, type := rep('Indexed cases', 7)]
setkey(index_age2, V1)
setnames(index_age2, c('age_group', 'n', 'type'))
index_age2[, age_group := factor(age_group, levels = c(1:7),
                                 labels = c('0-9',
                                            '10-19',
                                            '20-29',
                                            '30-39',
                                            '40-49',
                                            '50-59',
                                            '>=60'))]


in_degree <- sna::degree(net_before, cmode = 'indegree')
out_degree <- sna::degree(net_before, cmode = 'outdegree')
nonsingletons1 <- which(in_degree != 0 | out_degree != 0)
nonsingletons_age1 <- data.table(id = network::get.vertex.attribute(net_before, 'vertex.names')[nonsingletons1])
nonsingletons_age1[, age := plyr::mapvalues(id,
                                            from = COVID19_CJ_before$id1,
                                            to = COVID19_CJ_before$age_group1, warn_missing = F)]
nonsingletons_age1 <- as.data.table(table(nonsingletons_age1$age))

in_degree <- sna::degree(net_after, cmode = 'indegree')
out_degree <- sna::degree(net_after, cmode = 'outdegree')
nonsingletons2 <- which(in_degree != 0 | out_degree != 0)
nonsingletons_age2 <- data.table(id = network::get.vertex.attribute(net_after, 'vertex.names')[nonsingletons2])
nonsingletons_age2[, age := plyr::mapvalues(id,
                                            from = COVID19_CJ_after$id1,
                                            to = COVID19_CJ_after$age_group1, warn_missing = F)]
nonsingletons_age2 <- as.data.table(table(nonsingletons_age2$age))

index_age_period <- rbind(index_age1[, Period := rep('I', 7)], index_age2[, Period := rep('II', 7)])
index_age_period[Period == 'I', prop := n / nonsingletons_age1$N]
index_age_period[Period == 'II', prop := n / nonsingletons_age2$N]
generation_age_period[Period == 'I', prop := n / nonsingletons_age1$N]
generation_age_period[Period == 'II', prop := n / nonsingletons_age2$N]

age_index_last_period <- rbind(index_age_period, generation_age_period)

FigureS7a <- ggplot(age_index_last_period, aes(x = age_group, y = prop, fill = Period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    facet_wrap(~type) +
    geom_text(aes(label = scales::percent(prop, accuracy = 1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    # annotate('text', x = 1.75, y = 0.95, size = 10, family = 'Calibri', label = "p-value = 0.038",
    #          fontface = 'italic') +
    scale_fill_jama() +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1.05)) +
    labs(x = 'Age group', y = 'Proportion of infection') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(text = element_text(size = 30, family = FONTS),
          legend.position = 'top')

ggarrange(FigureS7a, font.label = list(size = 25, face = 'plain', family = FONTS), labels = '(a)',
          common.legend = TRUE, legend = "top")
ggsave(filename = 'figure/FigureS7a.svg', width = 16, height = 10)
ggsave(filename = 'figure/FigureS7a.pdf', width = 14, height = 7, device = cairo_pdf)

##########
in_degree <- sna::degree(net_before, cmode = 'indegree')
out_degree <- sna::degree(net_before, cmode = 'outdegree')
nonsingletons1 <- which(in_degree != 0 | out_degree != 0)
outdegree_age1 <- data.table(id = network::get.vertex.attribute(net_before, 'vertex.names')[nonsingletons1])
outdegree_age1[, age := plyr::mapvalues(id,
                                        from = COVID19_CJ_before$id1,
                                        to = COVID19_CJ_before$age_group1, warn_missing = F)]
outdegree_age1[, outdegree := out_degree[nonsingletons1]]

outdegree_age1 <- outdegree_age1 %>%
    group_by(age) %>%
    summarise(outdegree = mean(outdegree),
              n = n()) %>%
    as.data.table()

outdegree_age1[, age := factor(age, labels = c('0-9',
                                               '10-19',
                                               '20-29',
                                               '30-39',
                                               '40-49',
                                               '50-59',
                                               '>=60'))]

in_degree <- sna::degree(net_after, cmode = 'indegree')
out_degree <- sna::degree(net_after, cmode = 'outdegree')
nonsingletons2 <- which(in_degree != 0 | out_degree != 0)
outdegree_age2 <- data.table(id = network::get.vertex.attribute(net_after, 'vertex.names')[nonsingletons2])
outdegree_age2[, age := plyr::mapvalues(id,
                                        from = COVID19_CJ_after$id1,
                                        to = COVID19_CJ_after$age_group1, warn_missing = F)]
outdegree_age2[, outdegree := out_degree[nonsingletons2]]

outdegree_age2 <- outdegree_age2 %>%
    group_by(age) %>%
    summarise(outdegree = mean(outdegree),
              n = n()) %>%
    as.data.table()

outdegree_age2[, age := factor(age, labels = c('0-9',
                                               '10-19',
                                               '20-29',
                                               '30-39',
                                               '40-49',
                                               '50-59',
                                               '>=60'))]

outdegree_age_period <- rbind(outdegree_age1[, Period := rep('I', 7)],
                              outdegree_age2[, Period := rep('II', 7)])

FigureS7b <- ggplot(outdegree_age_period, aes(x = fct_inorder(age), y = outdegree, fill = Period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = round(outdegree, 2), y = outdegree),
              position = position_dodge(0.9),
              vjust = -0.5, size = 6, family = FONTS) +
    scale_fill_jama() +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2),
                       limits = c(0, 1.2)) +
    labs(x = 'Age group', y = 'Average out-degree') +
  guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'top', text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1))


#####################
household_dist_before <- trans_net3 %>%
    subset(id2 %in% plyr::mapvalues(trans_before_id$id2, from = COVID19$id2, to = 1:1349, warn_missing = F)) %>%
    group_by(is_last) %>%
    dplyr::summarise(ratio = mean(household),
                     n = n()) %>%
    as.data.table()

household_dist_after <- trans_net3 %>%
    subset(id2 %in% plyr::mapvalues(trans_after_id$id2, from = COVID19$id2, to = 1:1349, warn_missing = F)) %>%
    group_by(is_last) %>%
    dplyr::summarise(ratio = mean(household),
                     n = n()) %>%
    as.data.table()

household_dist_period <- rbind(household_dist_before[, Period := rep('I', 2)],
                               household_dist_after[, Period := rep('II', 2)])

household_dist_period[, is_last := factor(is_last, level = c(FALSE, TRUE),
                                          labels = c("FALSE", "TRUE"))]

FigureS7c <- ggplot(household_dist_period, aes(x = is_last, y = ratio, fill = Period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(ratio, accuracy = 0.1), y = ratio),
              position = position_dodge(0.9),
              vjust = -0.5, size = 6, family = FONTS) +
    scale_y_continuous(breaks = seq(0, 0.70, 0.2),
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 0.71)) +
    scale_fill_jama() +
    labs(x = 'Terminal transmission', y = 'Proportion within household') +
    guides(fill = guide_legend(title = 'Period')) +
    theme_bw() +
    theme(legend.position = 'top', text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1))

ggarrange(FigureS7b, FigureS7c, font.label = list(size = 25, face = 'plain', family = FONTS),
          labels = c('(b)', '(c)'), common.legend = TRUE, legend = "top", widths = c(10, 7),
          vjust=0.6)
ggsave(filename = 'figure/FigureS7bc.svg', width = 16, height = 10)
ggsave(filename = 'figure/FigureS7bc.pdf', width = 14, height = 7, device = cairo_pdf)

###


#########
expr <- p_value_expr(test_results[measures == "outdegree", p_value])
p1 <- ggplot(degree_dist, aes(x = fct_inorder(degree), y = prop, fill = period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    annotate('text', x = 3, y = 0.95, size = 10, label = expr,
             fontface = 'italic', family = FONTS) +
    #scale_fill_manual(values=col) +
    scale_fill_jama() +
    scale_y_continuous(breaks = seq(0, 1.00, 0.20),
                       labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1.00)) +
    labs(x = 'Out-degree', y = 'Proportion') +
    guides(fill = guide_legend(title = 'Period')) +
    # annotation_custom(ggplotGrob(p_degree_error), xmin = 2.7, xmax = 4.25, ymax = 0.90, ymin = 0.3) +
    theme_bw() +
    theme(text = element_text(size = 29, family = FONTS),
          legend.position = 'left',
          legend.title = element_text(size = 30),
          axis.title.y = element_blank(),
          panel.border = element_rect(size = 1.1))


expr <- p_value_expr(test_results[measures == "distance", p_value])
p2 <- ggplot(dist_dist, aes(x = fct_inorder(dist), y = prop, fill = period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    annotate('text', x = 3, y = 0.95, size = 10, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_fill_jama() +
    scale_y_continuous(breaks = seq(0, 1.00, 0.20),
                       labels = scales::percent_format(accuracy = 2),
                       limits = c(0, 1.00)) +
    labs(x = 'Distances between nodes', y = 'Proportion') +
    guides(fill = guide_legend(title = 'Period')) +
    # annotation_custom(ggplotGrob(p_dist_error),
    #                   xmin = 2.7, xmax = 4.25, ymax = 0.90, ymin = 0.3) +
    theme_bw() +
    theme(text = element_text(size = 29, family = FONTS),
          legend.position = 'left',
          legend.title = element_text(size = 30),
          axis.title.y = element_blank(),
          panel.border = element_rect(size = 1.1))

expr <- p_value_expr(test_results[measures == "betweenness", p_value])
p3 <- ggplot(betw_dist, aes(x = fct_inorder(betw), y = prop, fill = period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    annotate('text', x = 3, y = 0.95, size = 10, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_fill_jama() +
    scale_y_continuous(breaks = seq(0, 1.00, 0.20),
                       labels = scales::percent_format(accuracy = 2),
                       limits = c(0, 1.00)) +
    labs(x = 'Betweenness centrality', y = 'Proportion') +
    guides(fill = guide_legend(title = 'Period')) +
    # annotation_custom(ggplotGrob(p_size_error), xmin = 2.7, xmax = 4.25, ymax = 0.90, ymin = 0.3) +
    theme_bw() +
    theme(text = element_text(size = 29, family = FONTS),
          legend.position = 'left',
          legend.title = element_text(size = 30),
          axis.title.y = element_blank(),
          panel.border = element_rect(size = 1.1))

expr <- p_value_expr(test_results[measures == "diameter", p_value])
p4 <- ggplot(diam_dist, aes(x = fct_inorder(diam), y = prop, fill = period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    annotate('text', x = 3, y = 0.95, size = 10, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_fill_jama() +
    scale_y_continuous(breaks = seq(0, 1.00, 0.20),
                       labels = scales::percent_format(accuracy = 2),
                       limits = c(0, 1.00)) +
    labs(x = 'Diameter of clusters', y = 'Proportion') +
    guides(fill = guide_legend(title = 'Period')) +
    # annotation_custom(ggplotGrob(p_diam_error), xmin = 2.7, xmax = 4.25, ymax = 0.90, ymin = 0.3) +
    theme_bw() +
    theme(text = element_text(size = 29, family = FONTS),
          legend.position = 'left',
          legend.title = element_text(size = 30),
          axis.title.y = element_blank(),
          panel.border = element_rect(size = 1.1))

expr <- p_value_expr(test_results[measures == "cluster_size", p_value])
p5 <- ggplot(size_dist, aes(x = fct_inorder(size), y = prop, fill = period)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    annotate('text', x = 3, y = 0.95, size = 10, label = expr,
             fontface = 'italic', family = FONTS) +
    scale_fill_jama() +
    scale_y_continuous(breaks = seq(0, 1.00, 0.20),
                       labels = scales::percent_format(accuracy = 2),
                       limits = c(0, 1.00)) +
    labs(x = 'Size of clusters', y = 'Proportion') +
    guides(fill = guide_legend(title = 'Period')) +
    # annotation_custom(ggplotGrob(p_size_error), xmin = 2.7, xmax = 4.25, ymax = 0.90, ymin = 0.3) +
    theme_bw() +
    theme(text = element_text(size = 29, family = FONTS),
          legend.position = 'left',
          legend.title = element_text(size = 30),
          axis.title.y = element_blank(),
          panel.border = element_rect(size = 1.1))

figureS8 <- ggarrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3, common.legend = TRUE, legend = 'top',
                      labels = c('(a)', '(b)', '(c)', '(d)', '(e)'), vjust=0.6,
                      font.label = list(size = 25, face = 'plain', family = FONTS))

figureS8 <- annotate_figure(figureS8, 
                            left = text_grob("Proportion", face = 'plain', 
                                             size = 26, rot = 90, family = FONTS))
figureS8
ggsave(filename = 'figure/FigureS8.svg', width = 20, height = 12)
ggsave(filename = 'figure/FigureS8.pdf', width = 16, height = 16, device = cairo_pdf)


########import and removal
COVID19$remove_period <- as.Date(COVID19$confirm_date) - as.Date(COVID19$onset_date)
remove_period <- COVID19[, .(time_gap = mean(remove_period)), by = onset_date]
setkey(remove_period, 'onset_date')
remove_period$period <- c(rep('I', 16), rep('II', 31))

sp <- smooth.spline(remove_period$onset_date, remove_period$time_gap)

COVID19[, period := ifelse(onset_date <= "2020-01-23", "I", "II")]
period_table <- COVID19 %>%
    group_by(period) %>%
    summarise(mean = as.numeric(round(mean(remove_period), 2)),
              se = as.numeric(round(sd(remove_period) / sqrt(n()), 2))) %>%
    as.data.table()

period_table <- rbind(COVID19 %>% 
                          summarise(period="Overall",
                                    mean = as.numeric(round(mean(remove_period), 2)),
                                    se = as.numeric(round(sd(remove_period) / sqrt(n()), 2))) %>%
                          as.data.table(),
                      period_table)

pa <- ggplot() +
    geom_point(data = remove_period, aes(x = onset_date, y = time_gap, color = period), size = 3) +
    geom_line(aes(x = as.Date(sp$x, origin = '1970-01-01'), y = sp$y)) +
    geom_vline(xintercept = as.Date('2020-01-24'), linetype = "dashed",
               color = "black", size = 1.5) +
    annotate('text', x = as.Date('2020-01-24'),
             y = 25, size = 10, label = "Annoucement of the highest-level response", family = FONTS) +
    theme_bw() +
    labs(x = 'Date', y = 'Average removal period') +
    #scale_color_discrete("Period:", labels = c('I', "II")) +
    scale_color_jama() +
    scale_x_date(breaks = as.Date(c('2020-1-08', '2020-01-24', '2020-02-10', '2020-02-23')),
                 date_labels = "%b %d") +
    annotation_custom(gridExtra::tableGrob(period_table, rows = NULL, cols = c('Period', 'Avg. Removal Period', "Std. Err."),
                                           gridExtra::ttheme_default(base_size = 24,
                                                                     padding = unit(c(8, 5), 'mm'),
                                                                     base_family = FONTS),
                                           ),
                      xmin = as.Date('2020-01-27'), ymin = 12) +
    theme(text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1),
          legend.title = element_text(size = 30),
          legend.position = 'top') +
    guides(color = guide_legend(override.aes = list(size = 5), title = 'Period'))
pa

wuhan_date <- COVID19[, .(count = sum(wuhan, na.rm = TRUE)), by = onset_date]
setkey(wuhan_date, 'onset_date')

wuhan_date$period <- c(rep('I', 16), rep('II', 31))
# write.csv(wuhan_date, file = "data/wuhan_date.csv")
wuhan_date[, sum(count), by = period]

wuhan_date %>%
  group_by(period) %>%
  summarise(N = sum(count))

wuhan_table <- data.table(Period = c('Overall', 'I', 'II'), mean = c(239 + 206, 239, 206))
sp2 <- smooth.spline(wuhan_date$onset_date, wuhan_date$count)
pb <- ggplot() +
    geom_point(data = wuhan_date, aes(x = onset_date, y = count, color = period), size = 3) +
    geom_line(aes(x = as.Date(sp2$x, origin = '1970-01-01'), y = sp2$y)) +
    geom_vline(xintercept = as.Date('2020-01-24'), linetype = "dashed",
               color = "black", size = 1.5) +
    annotate('text', x = as.Date('2020-01-24'),
             y = 45, size = 10, label = "Annoucement of the highest-level response", family = FONTS) +
    theme_bw() +
    labs(x = 'Date', y = 'Number of imported cases') +
    #scale_color_discrete("Period:", labels = c('I', "II")) +
    scale_color_jama() +
    scale_x_date(breaks = as.Date(c('2020-1-08', '2020-01-24', '2020-02-10', '2020-02-23')),
                 date_labels = "%b %d") +
    annotation_custom(gridExtra::tableGrob(wuhan_table, rows = NULL, cols = c('Period', 'Number of import cases'),
                                           gridExtra::ttheme_default(base_size = 24,
                                                                     padding = unit(c(8, 5), 'mm'),
                                                                     base_family = FONTS)),
                      xmin = as.Date('2020-01-27'), ymin = 21) +
    theme(text = element_text(size = 30, family = FONTS),
          panel.border = element_rect(size = 1.1),
          legend.title = element_text(size = 30),
          legend.position = 'top') +
    guides(color = guide_legend(override.aes = list(size = 5), title = '时期'))
pb
ggarrange(pa, pb, nrow = 2, ncol=1, common.legend = TRUE, legend = 'top', labels = c('(a)', '(b)'),
          font.label = list(size = 25, face = 'plain', family = FONTS),
          vjust=0.6)
ggsave(filename = 'figure/FigureS9.svg', width = 18, height = 8)
ggsave(filename = 'figure/FigureS9.pdf', width = 14, height = 14, device = cairo_pdf)

############################COVID19 observed networks##############################
COVID19_CJ <- COVID19_CJ[infect == 1, c('id1', 'id2', 'onset_date2', 'infect', 'household')]
setnames(COVID19_CJ, c('id', 'contact', 'day', 'infect', 'household'))
COVID19_id_information <- data.table(id = 1:1349,
                                     infected_day = as.integer(as.Date(COVID19$onset_date) - as.Date('2020-01-05')))
COVID19_day_information <- data.table(day = 1:50,
                                      new_symp = sapply(as.Date('2020-01-05') + 0:49,
                                                        function(day, onset_date) {
                                                            return(sum(onset_date == day))
                                                        }, COVID19$onset_date))

observed_dynamics <- rbindlist(
    lapply(4:50, function(i) { return(sum_graph_dt_period_covid(i, 0, COVID19_CJ, COVID19_day_information,
                                                                COVID19_id_information, 1349)) }))

observed_dynamics <- melt(observed_dynamics, id.vars = 'period', variable.name = 'Measures')

observed_dynamics <- rbind(observed_dynamics,
                           data.table(period = 4:50, Measures = "p_household",
                                      value = sapply(4:50, function(i, x) { return(mean(x[day <= i]$household)) },
                                                     COVID19_CJ)))
observed_dynamics[, Measures := factor(Measures,
                                       levels = unique(observed_dynamics$Measures),
                                       labels = c('A: Number of the new-onset',
                                                  'B: Number of singletons',
                                                  'C: Percentage of singletons',
                                                  'D: Average outdegree',
                                                  'E: Average diameter of clusters',
                                                  'F: Average size of clusters',
                                                  'G: Average shortest path length',
                                                  'H: Proportion within household'))]
setnames(observed_dynamics, c('day', 'Measures', 'median'))
observed_dynamics[, min := NA]
observed_dynamics[, max := NA]
observed_dynamics[, scenario := 'Observed']
setcolorder(observed_dynamics, c('day', 'Measures', 'min', 'median', 'max', 'scenario'))
observed_dynamics[, day := day - 3]

scenario1_all_contact_infect <- data.table(read_parquet('result1/all_contact_infect.parquet'))[scenario == 1]
scenario1_day_information <- data.table(read_parquet('result1/day_information.parquet'))[scenario == 1]
scenario1_id_information <- data.table(read_parquet('result1/id_information.parquet'))[scenario == 1]

Rt <- Rt_dt(1, 200, scenario1_day_information, scenario1_all_contact_infect,
            scenario1_id_information, wuhan_date)
simulation_summary_graph_total <- Rt

simulation_summary_graph_python <- fread('result1/simulation_summary_total.csv')[scenario == 1]
simulation_summary_graph_total <- rbind(simulation_summary_graph_python, simulation_summary_graph_total)
# simulation_summary_graph_total[scenario == 1, scenario := as.factor(scenario)]

write.csv(simulation_summary_graph_total, file = 'result2/simulation_summary_total1.csv')

scenario1_summary <- fread(file = "result2/simulation_summary_total1.csv")[scenario == 1]
scenario1_summary[, V1 := NULL]
scenario1_summary[, scenario := rep('Reconstructed', nrow(scenario1_summary))]
scenario1_summary[Measures == 'Percentage of infection', Measures := 'Number of infection']
scenario1_summary[Measures == 'Number of infection', c('min', 'median', 'max')] <-
    scenario1_summary[Measures == 'Number of infection', lapply(.SD, function(i) i * 20000),
                      .SDcols = c('min', 'median', 'max')]
scenario1_summary <- scenario1_summary[!Measures %in% c('Number of infection', 'Effective reproduction numbers'),]
# The first day is day 8.
scenario1_summary <- scenario1_summary[day <= 47,]
scenario1_summary[, Measures := factor(Measures,
                                       levels = unique(scenario1_summary$Measures),
                                       labels = c('A: Number of the new-onset',
                                                  'B: Proportion within household',
                                                  'C: Average outdegree',
                                                  'D: Average shortest path length',
                                                  'E: Average betweenness',
                                                  'F: Average diameter of clusters',
                                                  'G: Average size of clusters',
                                                  'H: Number of singletons',
                                                  'I: Percentage of singletons'))]

obs_vs_recon <- rbind(observed_dynamics, scenario1_summary)
# simulate_dynamics <- ggplot(data = obs_vs_recon, aes(color = scenario, linetype = scenario)) +
#     geom_line(aes(x = day, y = median), size = 0.9) +
#     geom_ribbon(aes(x = day, ymin = min, ymax = max, fill = scenario), alpha = .3) +
#     facet_wrap_custom(~Measures, ncol = 2, scale = 'free_y',
#                       scale_overrides = list(scale_override(3,
#                                                             scale_y_continuous(labels = scales::percent_format())),
#                                              scale_override(8,
#                                                             scale_y_continuous(labels = scales::percent_format())))) +
#     scale_fill_manual(values = c("#d95f0e", "#2c7fb8")) +
#     theme_bw() +
#     labs(color = 'Scenario', linetype = 'Scenario', fill = 'Scenario') +
#     theme(text = element_text(size = 25),
#           legend.position = 'top',
#           legend.title = element_blank(),
#           axis.title.y = element_blank(),
#           panel.border = element_rect(size = 1.1)) +
#     guides(color = guide_legend(override.aes = list(size = 0.7),
#                                 nrow = 2, byrow = T))
obs_vs_recon[, scenario := factor(scenario, levels = c("Observed", "Reconstructed"),
                                  labels = c("Observed", "Re-simulated"))]
simulate_dynamics <- ggplot(data = obs_vs_recon[Measures == "A: Number of the new-onset"], 
                            aes(color = scenario, linetype = scenario)) +
    geom_line(aes(x = day, y = median), size = 0.9) +
    geom_ribbon(aes(x = day, ymin = min, ymax = max, fill = scenario), alpha = .3) +
    scale_fill_manual(values = c("#d95f0e", "#2c7fb8")) +
    theme_bw() +
    labs(color = 'Scenario', linetype = 'Scenario', fill = 'Scenario', 
         x = "Day", y = "Number of the new-onset") +
    theme(text = element_text(size = 25, family = FONTS),
          legend.position = 'top',
          legend.title = element_blank(),
          panel.border = element_rect(size = 1.1)) +
    guides(color = guide_legend(override.aes = list(size = 0.7),
                                nrow = 2, byrow = T))

simulate_dynamics
ggsave(filename = "figure/FigureS10.svg", width = 16, height = 8)
ggsave(filename = "figure1/FigureS10.pdf", width = 16, height = 8, device = cairo_pdf)

############################COVID19 curves##############################
curves <- observed_dynamics[Measures == 'A: Number of the new-onset']
curves[, day := seq(0, nrow(curves) - 1) + as.Date('2020-01-08')]
ggplot(data = curves) + 
  geom_bar(aes(x = day, y = median), stat='identity') + 
  annotate("rect", xmin = as.Date('2020-01-08') + 16,
           xmax = as.Date('2020-01-08') + 33, ymin = 0, ymax = Inf, alpha = .2) +
  annotate('text', x = as.Date('2020-01-07') + 25,
           y = 85, size = 7, label = "一级管控期间", family = FONTS) +
  theme_bw() +
  scale_x_date(breaks = as.Date(c('2020-1-08', '2020-01-24', '2020-02-10', '2020-02-23')),
               date_labels = "%b %d") +
  scale_y_continuous(breaks = seq(0, 100, 25), limit=c(0, 100)) +
  labs(x = "日期", y = "每日新增有症状病例数") +
  theme(text = element_text(size = 25, family = FONTS),
        legend.title = element_blank(),
        panel.border = element_rect(size = 1.1))

ggsave(filename = "figure/curves.pdf", width = 16, height = 8, device = cairo_pdf)



############################age distribution##############################
age_census <- fread("data/ZJ_pop.csv")
age_census[, age := fct_inorder(c('0-4', '5-9', '10-14', '15-19',
                                  '20-24', '25-29', '30-34', '35-39',
                                  '40-44', '45-49', '50-54', '55-59',
                                  '60-64', '>=65'))][, prop := freq / sum(freq)]


ggplot(age_census, aes(x = fct_inorder(age), y = prop)) +
    geom_bar(stat = 'identity', position = position_dodge(0.9)) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1), y = prop),
              position = position_dodge(0.9), vjust = -0.5, size = 6, family = FONTS) +
    scale_y_continuous(breaks = seq(0, 0.16, 0.04),
                       labels = scales::percent_format(accuracy = 2),
                       limits = c(0, 0.16)) +
    labs(x = 'Age group', y = 'Proportion') +
    theme_bw() +
    theme(text = element_text(size = 30, family = FONTS),
          axis.text.x = element_text(angle = 45, vjust = 0.6, hjust=0.5, 
                                     size = 24, family = FONTS),
          legend.position = 'left',
          legend.title = element_text(size = 30),
          panel.border = element_rect(size = 1.1))

ggsave(filename = "figure/age_distribution.pdf", width = 15, height = 6, device = cairo_pdf)
