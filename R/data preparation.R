library(ggplot2)
library(ggsci)
library(forcats)
library(data.table)
library(extrafont)
library(ggpubr)
library(viridis)
library(GGally)
library(network)
library(dplyr)
library(latex2exp)
library(arrow)

##
COVID19 <- data.table(readxl::read_excel('data/COVID19.xlsx'))
COVID19$onset_date <- as.Date(COVID19$onset_date) + 1
COVID19$confirm_date <- as.Date(COVID19$confirm_date) + 1
COVID19$import_date <- as.Date(COVID19$import_date)
COVID19 <- COVID19[onset_date > '2020-01-05']

wuhan_date <- COVID19 %>%
  group_by(onset_date) %>%
  summarise(count=sum(wuhan, na.rm = TRUE)) %>% as.data.table()
# 
write.csv(wuhan_date, "data/wuhan_date.csv")

wuhan_date <- fread('data/wuhan_date.csv')
wuhan_date[, V1 := NULL]
wuhan_date[, cumcount := cumsum(count)]

import_cases <- COVID19[wuhan == 1, c("onset_date", "import_date")]
import_cases$onset_date <- import_cases$onset_date - as.Date('2020-01-08')
import_cases$import_date <- import_cases$import_date - as.Date('2020-01-08')
write.csv(import_cases, "data/import_cases.csv")

######################COVID19 observed network preparation##########################
set.seed(2021)
COVID19$parent <- apply(COVID19, MARGIN = 1, FUN = sample_parent)
COVID19_CJ <- CJ(id1 = COVID19$id2, id2 = COVID19$id2)

COVID19_CJ <- COVID19_CJ %>%
  merge(COVID19[, .(id2 = id2, family1 = family_index2, onset_date1 = onset_date)],
        by.x = 'id1', by.y = 'id2', sort = FALSE) %>%
  merge(COVID19[, .(id2 = id2, family2 = family_index2, onset_date2 = onset_date, parent = parent)],
        by = 'id2', sort = FALSE)

COVID19_CJ[, onset_date1 := as.integer(onset_date1 - as.Date('2020-01-04'))]
COVID19_CJ[, onset_date2 := as.integer(onset_date2 - as.Date('2020-01-04'))]

setcolorder(COVID19_CJ, c('id1', 'id2', 'family1', 'family2', 'onset_date1', 'onset_date2', 'parent'))
COVID19_CJ$infect <- ifelse(!is.na(COVID19_CJ$parent) & (COVID19_CJ$id1 == COVID19_CJ$parent), 1, 0)
COVID19_CJ[, household := ifelse(!is.na(family1) &
                                   !is.na(family2) &
                                   family1 == family2, 1, 0)]
COVID19_CJ[, parent := NULL]
# system.time(COVID19_CJ[, age_group1 := plyr::mapvalues(id1, from = COVID19$id2, to = COVID19$age_year, warn_missing = F)])
# COVID19_CJ[, age_group2 := plyr::mapvalues(id2, from = COVID19$id2, to = COVID19$age_year, warn_missing = F)]
# COVID19_CJ[, family1 := plyr::mapvalues(id1, from = COVID19$id2, to = COVID19$family_index2, warn_missing = F)]
# COVID19_CJ[, family2 := plyr::mapvalues(id2, from = COVID19$id2, to = COVID19$family_index2, warn_missing = F)]

COVID19_CJ <- COVID19_CJ %>%
  merge(COVID19[, .(id2 = id2, age_group1 = age_year)], by.x = "id1", by.y = "id2", sort = FALSE) %>%
  merge(COVID19[, .(id2 = id2, age_group2 = age_year)], by.x = "id2", by.y = "id2", sort = FALSE)

setcolorder(COVID19_CJ, c('id1', 'id2', 'onset_date1', 'onset_date2', 'infect', 'household',
                          'age_group1', 'age_group2', 'family1', 'family2'))

COVID19_CJ[, age_group1 := cut(as.numeric(age_group1), breaks = c(seq(0, 60, 10), 100), right = FALSE, labels = FALSE)]
COVID19_CJ[, age_group2 := cut(as.numeric(age_group2), breaks = c(seq(0, 60, 10), 100), right = FALSE, labels = FALSE)]

COVID19_CJ[, id1 := COVID19_CJ[, "id1"] %>%
             merge(COVID19[, .(id2 = id2, num = 1:1349)], by.x = "id1", by.y = "id2", sort = FALSE) %>%
             subset(select = "num")]
COVID19_CJ[, id2 := COVID19_CJ[, "id2"] %>%
             merge(COVID19[, .(id2 = id2, num = 1:1349)], by = "id2", sort = FALSE) %>%
             subset(select = "num")]


###################COVID19 observed network construction#######################
trans_net <- dcast(COVID19_CJ, id1 ~ id2, value.var = 'infect')
adj_net <- as.matrix(trans_net[, -1])
net <- network(adj_net)
network::set.vertex.attribute(net, 'age_interval',
                              cut(as.numeric(plyr::mapvalues(colnames(adj_net),
                                                             from = 1:1349, to = COVID19$age_year)),
                                  breaks = c(seq(0, 60, 10), 100),
                                  right = FALSE,
                                  labels = FALSE))

household_wide <- as.matrix((dcast(COVID19_CJ, id1 ~ id2, value.var = 'household'))[, -1])
network::set.edge.value(net, 'household', ifelse(household_wide > 0, "#e41a1c", "grey75"))


origin <- COVID19[, .(origin_time = min(onset_date)), by = cluster]
# COVID19[, origin_time :=
#               as.Date(plyr::mapvalues(cluster, from = origin$cluster,
#                                       to = origin$origin_time, warn_missing = FALSE), origin = '1970-01-01')]
# COVID19_CJ[, origin_time := as.Date(plyr::mapvalues(id1, from = 1:1349, to = COVID19$origin_time), origin = '1970-01-01')]

COVID19 <- COVID19 %>%
  merge(origin, by = "cluster", sort = FALSE)

COVID19[is.na(cluster), origin_time := onset_date]

COVID19_CJ <- COVID19_CJ %>%
  merge(COVID19[, .(id1 = 1:1349, origin_time = origin_time)], by = "id1", sort = FALSE)

T_graph <- intergraph::asIgraph(net)
dist_T <- igraph::distances(T_graph, mode = 'out')
colnames(dist_T) <- rownames(dist_T) <- 1:1349

in_degree <- sna::degree(net, cmode = 'indegree')
out_degree <- sna::degree(net, cmode = 'outdegree')
indexed <- which(in_degree == 0 & out_degree > 0)

a <- dist_T[indexed,]

distance_dt <- data.table(reshape2::melt(a))
setnames(distance_dt, c('id1', 'id2', 'distance'))
distance_dt <- distance_dt[distance != Inf]
distance_dt <- distance_dt[distance != 0]

a <- distance_dt[, list(max_generation = max(distance)), by = id1]
distance_dt[, max_generation := plyr::mapvalues(id1, from = a$id1, to = a$max_generation)]


COVID19_CJ[, generation := 0]
a <- as.integer(plyr::mapvalues(COVID19_CJ[id2 %in% distance_dt$id2, id2], from = distance_dt$id2,
                                to = distance_dt$distance))
COVID19_CJ[id2 %in% distance_dt$id2, generation := a]

COVID19_before <- COVID19[origin_time <= TIME_CUT]
COVID19_CJ_before <- CJ(id1 = COVID19_before$id2, id2 = COVID19_before$id2)
COVID19_CJ_before <- COVID19_CJ_before %>%
  merge(COVID19_before[, .(id2 = id2, family1 = family_index2, onset_date1 = onset_date)],
        by.x = 'id1', by.y = 'id2', sort = FALSE) %>%
  merge(COVID19_before[, .(id2 = id2, family2 = family_index2, onset_date2 = onset_date, parent = parent)],
        by = 'id2', sort = FALSE)

setcolorder(COVID19_CJ_before, c('id1', 'id2', 'family1', 'family2', 'onset_date1', 'onset_date2', 'parent'))
COVID19_CJ_before[, infect := ifelse(!is.na(parent) & (id1 == parent), 1, 0)]
COVID19_CJ_before[, household := ifelse(!is.na(family1) &
                                          !is.na(family2) &
                                          (family1 == family2), 1, 0)]
COVID19_CJ_before[, c('onset_date1', 'parent') := NULL]

COVID19_CJ_before <- COVID19_CJ_before %>%
  merge(COVID19_before[, .(id2 = id2, age_group1 = age_year)], by.x = "id1", by.y = "id2", sort = FALSE) %>%
  merge(COVID19_before[, .(id2 = id2, age_group2 = age_year)], by.x = "id2", by.y = "id2", sort = FALSE)

COVID19_CJ_before[, age_group1 := cut(as.numeric(age_group1), breaks = c(seq(0, 60, 10), 100),
                                      right = FALSE, labels = FALSE)]
COVID19_CJ_before[, age_group2 := cut(as.numeric(age_group2), breaks = c(seq(0, 60, 10), 100),
                                      right = FALSE, labels = FALSE)]

COVID19_CJ_before[, onset_date2 := as.integer(onset_date2 - as.Date('2020-01-04'))]

setcolorder(COVID19_CJ_before, c('id1', 'id2', 'onset_date2', 'infect', 'household',
                                 'age_group1', 'age_group2', 'family1', 'family2'))

trans_net_matrix_before <- dcast(COVID19_CJ_before, id1 ~ id2, value.var = 'infect')
adj_net <- as.matrix(trans_net_matrix_before[, -1])
net_before <- network(adj_net)

network::set.vertex.attribute(net_before, 'age_interval',
                              cut(as.numeric(plyr::mapvalues(colnames(adj_net),
                                                             from = COVID19_before$id2, to = COVID19_before$age_year)),
                                  breaks = c(seq(0, 60, 10), 100),
                                  right = FALSE,
                                  labels = FALSE))

household_wide <- as.matrix((dcast(COVID19_CJ_before, id1 ~ id2, value.var = 'household'))[, -1])
network::set.edge.value(net_before, 'household', ifelse(household_wide > 0, "#e41a1c", "grey75"))


COVID19_after <- COVID19[origin_time > TIME_CUT]
COVID19_CJ_after <- CJ(id1 = COVID19_after$id2, id2 = COVID19_after$id2)

COVID19_CJ_after <- COVID19_CJ_after %>%
  merge(COVID19_after[, .(id2 = id2, family1 = family_index2, onset_date1 = onset_date)],
        by.x = 'id1', by.y = 'id2', sort = FALSE) %>%
  merge(COVID19_after[, .(id2 = id2, family2 = family_index2, onset_date2 = onset_date, parent = parent)],
        by = 'id2', sort = FALSE)

setcolorder(COVID19_CJ_after, c('id1', 'id2', 'family1', 'family2', 'onset_date1', 'onset_date2', 'parent'))
COVID19_CJ_after[, infect := ifelse(!is.na(parent) & (id1 == parent), 1, 0)]
COVID19_CJ_after[, household := ifelse(!is.na(family1) &
                                         !is.na(family2) &
                                         (family1 == family2), 1, 0)]
COVID19_CJ_after[, c('onset_date1', 'parent') := NULL]

COVID19_CJ_after <- COVID19_CJ_after %>%
  merge(COVID19_after[, .(id2 = id2, age_group1 = age_year)], by.x = "id1", by.y = "id2", sort = FALSE) %>%
  merge(COVID19_after[, .(id2 = id2, age_group2 = age_year)], by.x = "id2", by.y = "id2", sort = FALSE)

COVID19_CJ_after[, age_group1 := cut(as.numeric(age_group1), breaks = c(seq(0, 60, 10), 100),
                                     right = FALSE, labels = FALSE)]
COVID19_CJ_after[, age_group2 := cut(as.numeric(age_group2), breaks = c(seq(0, 60, 10), 100),
                                     right = FALSE, labels = FALSE)]
COVID19_CJ_after[, onset_date2 := as.integer(onset_date2 - as.Date('2020-01-04'))]

setcolorder(COVID19_CJ_after, c('id1', 'id2', 'onset_date2', 'infect', 'household',
                                'age_group1', 'age_group2', 'family1', 'family2'))

trans_net_matrix_after <- dcast(COVID19_CJ_after, id1 ~ id2, value.var = 'infect')
adj_net <- as.matrix(trans_net_matrix_after[, -1])
net_after <- network(adj_net)

network::set.vertex.attribute(net_after, 'age_interval',
                              cut(as.numeric(plyr::mapvalues(colnames(adj_net),
                                                             from = COVID19_after$id2, to = COVID19_after$age_year)),
                                  breaks = c(seq(0, 60, 10), 100),
                                  right = FALSE,
                                  labels = FALSE))

household_wide <- as.matrix((dcast(COVID19_CJ_after, id1 ~ id2, value.var = 'household'))[, -1])
network::set.edge.value(net_after, 'household', ifelse(household_wide > 0, "#e41a1c", "grey75"))