TIME_CUT <- '2020-01-20'

source('R/utils.R')
source('R/facet_wrap_custom.R')
source("R/data preparation.R")
source("R/statistical test.R")
test_results

TIME_CUT <- '2020-02-01'

source('R/utils.R')
source('R/facet_wrap_custom.R')
source("R/data preparation.R")
source("R/statistical test.R")
test_results


sensitivity_parent <- function(seed){
  test_to_df <- function(measures, t, ci1, ci2){
    if (is.null(ci1)){
      data.table(measures=measures,
                 period1=NaN,
                 period2=NaN,
                 period1_ci=NaN,
                 period2_ci=NaN,
                 p_value=as.numeric(t$p.value))
    }else{
      data.table(measures=measures,
                 period1=t$estimate[1],
                 period2=t$estimate[2],
                 period1_ci=paste("(", paste(round(ci1$percent[4:5], 3), collapse=", "), ")", sep=""),
                 period2_ci=paste("(", paste(round(ci2$percent[4:5], 3), collapse=", "), ")", sep=""),
                 p_value=as.numeric(t$p.value))
    }
  }
  
  boot.ci2 <- function(boot_samples, type="perc"){
    if (type=="perc"){
      name <- "percent"
    }else if(type == 'norm'){
      name <- "normal"
    }else{
      name <- "CI"
    }
    if (is.null(boot.ci(boot_samples, type=type))){
      tmp <- list(to=boot_samples$t0, name=c(0.95, 250, 9750, boot_samples$t0, boot_samples$t0))
      names(tmp)[2] <- name
      return(tmp)
    }else{
      return(boot.ci(boot_samples, type=type))
    }
  }
  N_TIMES <- 100
  test_results <- data.table()
  set.seed(seed)
  COVID19 <- data.table(readxl::read_excel('data/COVID19.xlsx'))
  COVID19$onset_date <- as.Date(COVID19$onset_date) + 1
  COVID19$confirm_date <- as.Date(COVID19$confirm_date) + 1
  COVID19$import_date <- as.Date(COVID19$import_date)
  COVID19 <- COVID19[onset_date > '2020-01-05']
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
  
  ####################test for proportion of singletons####################
  in_degree <- sna::degree(net, cmode = 'indegree')
  out_degree <- sna::degree(net, cmode = 'outdegree')
  nonsingletons <- which(in_degree != 0 | out_degree != 0)
  is_singletons <- data.table(id = 1:1349)
  is_singletons[, is_single := TRUE]
  is_singletons[id %in% nonsingletons, is_single := FALSE]
  
  
  is_singletons[, period := 2]
  is_singletons[id %in% (1:1349)[COVID19$origin_time <= TIME_CUT], period := 1]
  is_singletons[, out_degree := out_degree]
  
  singletons <- is_singletons %>%
    group_by(period) %>%
    summarise(is_singletons = sum(is_single),
              n_cases = n()) %>%
    as.data.table()
  
  # CI for proportion
  period1_ci <- boot.ci(boot(is_singletons[period == 1, is_single], bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci(boot(is_singletons[period == 2, is_single], bootMean, R = N_TIMES), type = "perc")
  
  test <- prop.test(singletons$is_singletons, singletons$n_cases, alternative = 'less')
  
  test_results <- rbind(test_results, test_to_df("prop_singletons", test, period1_ci, period2_ci))
  
  ####################test for proportion of super-spreaders####################
  super_spreaders <- is_singletons %>%
    subset(is_single == FALSE) %>%
    group_by(period) %>%
    summarise(super = sum(out_degree >= 3),
              n_cases = n(),
              prop = super / n_cases) %>%
    as.data.table()
  
  # CI for proportion
  period1_ci <- boot.ci(boot(is_singletons[is_single == FALSE & period == 1, out_degree >= 3], bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci(boot(is_singletons[is_single == FALSE & period == 2, out_degree >= 3], bootMean, R = N_TIMES), type = "perc")
  
  test <- prop.test(super_spreaders$super, super_spreaders$n_cases, alternative = 'greater')
  
  test_results <- rbind(test_results, test_to_df("prop_superspreaders", test, period1_ci, period2_ci))
  
  
  ####################test for out-degree across periods####################
  in_degree <- sna::degree(net, cmode = 'indegree')
  out_degree <- sna::degree(net, cmode = 'outdegree')
  nonsingletons <- which(in_degree != 0 | out_degree != 0)
  degree_period <- data.table(id = 1:1349, out_degree = out_degree)
  degree_period[, period := 2]
  degree_period[id %in% (1:1349)[COVID19$origin_time <= TIME_CUT], period := 1]
  degree_period <- degree_period[id %in% nonsingletons]
  
  # CI for out-degree by period
  degree_period %>%
    group_by(period) %>%
    summarise(mean_outdegree = mean(out_degree))
  period1_ci <- boot.ci(boot(degree_period[period == 1, out_degree], bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci(boot(degree_period[period == 2, out_degree], bootMean, R = N_TIMES), type = "perc")
  
  test <- t.test(out_degree ~ period, data = degree_period, alternative = 'greater')
  
  test_results <- rbind(test_results, test_to_df("avg_outdegree", test, period1_ci, period2_ci))
  
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
  
  test <- chisq.test(matrix(degree_dist$N, nrow = 2))
  
  test_results <- rbind(test_results, test_to_df("outdegree", test, NULL, NULL))
  
  ####################test for diameter across periods####################
  
  # CI for diameter by period
  period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = 'diameter'), bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci2(boot(summary_graph_feature(net_after, feature = 'diameter'), bootMean, R = N_TIMES), type = "perc")
  
  #
  test <- t.test(x = summary_graph_feature(net_before, feature = 'diameter'),
                 y = summary_graph_feature(net_after, feature = 'diameter'),
                 alternative = 'greater')
  
  test_results <- rbind(test_results, test_to_df("avg_diameter", test, period1_ci, period2_ci))
  
  diam_period <- rbind(data.table(period = 1, dist = summary_graph_feature(net_before, feature = 'diameter')),
                       data.table(period = 2, dist = summary_graph_feature(net_after, feature = 'diameter')))
  
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
  
  test <- chisq.test(matrix(diam_dist$N, nrow = 2))
  
  test_results <- rbind(test_results, test_to_df("diameter", test, NULL, NULL))
  
  ####################test for distance across periods####################
  # CI for distance by period
  period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = 'distance'), bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci2(boot(summary_graph_feature(net_after, feature = 'distance'), bootMean, R = N_TIMES), type = "perc")
  
  test <- t.test(x = summary_graph_feature(net_before, feature = 'distance'),
                 y = summary_graph_feature(net_after, feature = 'distance'),
                 alternative = 'greater')
  
  test_results <- rbind(test_results, test_to_df("avg_distance", test, period1_ci, period2_ci))
  
  
  data.table(diameter = summary_graph_feature(net_before, feature = 'distance'), period = 1) %>%
    rbind(data.table(diameter = summary_graph_feature(net_after, feature = 'distance'), period = 2)) %>%
    group_by(period) %>%
    table() %>%
    # chisq.test()
    fisher.test()
  
  
  dist_period <- rbind(data.table(period = 1, dist = summary_graph_feature(net_before, feature = 'distance')),
                       data.table(period = 2, dist = summary_graph_feature(net_after, feature = 'distance')))
  
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
  
  
  test <- chisq.test(matrix(dist_dist$N, nrow = 2))
  
  test_results <- rbind(test_results, test_to_df("distance", test, NULL, NULL))
  
  ####################test for cluster sizes across periods####################
  
  # CI for cluster sizes by period
  period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = 'cluster size'), bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci(boot(summary_graph_feature(net_after, feature = 'cluster size'), bootMean, R = N_TIMES), type = "perc")
  
  test <- t.test(x = summary_graph_feature(net_before, feature = 'cluster size'),
                 y = summary_graph_feature(net_after, feature = 'cluster size'),
                 alternative = 'greater')
  
  test_results <- rbind(test_results, test_to_df("avg_cluster_size", test, period1_ci, period2_ci))
  
  size_period <- rbind(data.table(period = 1, dist = summary_graph_feature(net_before, feature = 'cluster size')),
                       data.table(period = 2, dist = summary_graph_feature(net_after, feature = 'cluster size')))
  
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
  
  
  test <- chisq.test(matrix(size_dist$N, nrow = 2))
  
  test_results <- rbind(test_results, test_to_df("cluster size", test, NULL, NULL))
  
  ####################test for betweenness across periods####################
  # overall interval
  feature <- "betweenness"
  
  # CI for cluster sizes by period
  period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = feature), bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci(boot(summary_graph_feature(net_after, feature = feature), bootMean, R = N_TIMES), type = "perc")
  
  test <- t.test(x = summary_graph_feature(net_before, feature = feature),
                 y = summary_graph_feature(net_after, feature = feature))
  
  test_results <- rbind(test_results, test_to_df("avg_betweenness", test, period1_ci, period2_ci))
  
  betw_period <- rbind(data.table(period = 1, betw = summary_graph_feature(net_before, feature = feature)),
                       data.table(period = 2, betw = summary_graph_feature(net_after, feature = feature)))
  
  betw_dist <- betw_period %>%
    group_by(period) %>%
    summarise('0' = sum(betw == 0),
              '1' = sum(betw == 1),
              '2' = sum(betw == 2),
              '>=3' = sum(betw >= 3)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "betw") %>%
    as.data.table()
  betw_dist[, total := sum(N), by = period][, prop := N / total][, total := NULL]
  betw_dist[, period := factor(period, levels = c(1, 2), labels = c("I", "II"))]
  
  
  test <- chisq.test(matrix(betw_dist$N, nrow = 2))
  
  test_results <- rbind(test_results, test_to_df("betweenness", test, NULL, NULL))
  
  
  ####################test for household proportion across periods####################
  household_before <- COVID19_CJ_before[infect != 0, c('id1', 'id2', 'household')]
  household_before[is.na(household), household := 0]
  household_after <- COVID19_CJ_after[infect != 0, c('id1', 'id2', 'household')]
  household_after[is.na(household), household := 0]
  
  # overall interval for household
  household <- COVID19_CJ[infect != 0, ]$household
  household[is.na(household)] <- 0
  
  # CI for household proportion by period
  period1_ci <- boot.ci(boot(household_before$household, bootMean, R = N_TIMES), type = "perc")
  period2_ci <- boot.ci(boot(household_after$household, bootMean, R = N_TIMES), type = "perc")
  
  test <- prop.test(x = c(sum(household_before$household), sum(household_after$household)),
                    n = c(nrow(household_before), nrow(household_after)), alternative = "less")
  
  test_results <- rbind(test_results, test_to_df("prop_household", test, period1_ci, period2_ci))
  
  ####################test for removal periods####################
  COVID19$remove_period <- as.Date(COVID19$confirm_date) - as.Date(COVID19$onset_date)
  remove_period <- COVID19[, .(time_gap = mean(remove_period)), by = onset_date]
  setkey(remove_period, 'onset_date')
  remove_period$period <- c(rep('I', 16), rep('II', 31))
  
  # write.csv(remove_period, file="data/remove_period.csv")
  
  test <- notrend_test(as.numeric(remove_period$time_gap), test = "MK")
  
  test_results <- rbind(test_results, test_to_df("removal_trend", test, NULL, NULL))
  
  remove_period_period <- data.table(id = 1:1349, remove_period = COVID19$remove_period,
                                     is_symp = ifelse(COVID19$case_dignosis2 == 'confirmed_cases', 1, 0))
  remove_period_period[, period := 2]
  remove_period_period[id %in% (1:1349)[COVID19$onset_date <= TIME_CUT], period := 1]
  
  test <- t.test(data = remove_period_period, remove_period ~ period, alternative = 'greater')
  
  test_results <- rbind(test_results, test_to_df("avg_removal", test, NULL, NULL))
  
  COVID19 %>%
    subset(onset_date >= '2020-02-10') %>%
    summarise(mean_remove=mean(remove_period),
              se_remove=sd(remove_period) / sqrt(n()))
  
  ####################multiple adjustment####################
  test_results[, p_value := round(p.adjust(p_value, method = "BH"), 4)]
  return(test_results)
}

parent_results <- lapply(1:100, sensitivity_parent)

a <- rowSums(sapply(parent_results, function(x) {(x$p_value <= 0.05)})) / 100

parent_results[[1]]$measures[which(! a %in% c(0, 1))]

