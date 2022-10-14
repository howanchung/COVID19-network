N_TIMES <- 100
USE_BOOT <- TRUE

test_results <- data.table()

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

# num of clusters
cluster <- igraph::decompose(intergraph::asIgraph(net))
cluster <- cluster[sapply(cluster, igraph::vcount) > 1]
cluster_before <- igraph::decompose(intergraph::asIgraph(net_before))
cluster_before <- cluster_before[sapply(cluster_before, igraph::vcount) > 1]
cluster_after <- igraph::decompose(intergraph::asIgraph(net_after))
cluster_after <- cluster_after[sapply(cluster_after, igraph::vcount) > 1]
print(length(cluster))
print(length(cluster_before))
print(length(cluster_after))


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

# overall interval
mean(is_singletons$is_single)
set.seed(2021)
boot.ci(boot(is_singletons$is_single, bootMean, R = N_TIMES), type = "perc")

singletons <- is_singletons %>%
  group_by(period) %>%
  summarise(is_singletons = sum(is_single),
            n_cases = n()) %>%
  as.data.table()

singletons

# CI for proportion
period1_ci <- boot.ci(boot(is_singletons[period == 1, is_single], bootMean, R = N_TIMES), type = "perc")
period2_ci <- boot.ci(boot(is_singletons[period == 2, is_single], bootMean, R = N_TIMES), type = "perc")

test <- prop.test(singletons$is_singletons, singletons$n_cases, alternative = 'less')

test_results <- rbind(test_results, test_to_df("prop_singletons", test, period1_ci, period2_ci))



# CI for measures
# CI for proportion
get_number_super_spreaders <- function(net_obj){
  in_degree <- sna::degree(net_obj, cmode = 'indegree')
  out_degree <- sna::degree(net_obj, cmode = 'outdegree')
  nonsingletons <- which(in_degree != 0 | out_degree != 0)
  return(c(sum(out_degree[nonsingletons] >= 3), length(nonsingletons)))
}
prop_super_spreaders <- function(net_obj){
  in_degree <- sna::degree(net_obj, cmode = 'indegree')
  out_degree <- sna::degree(net_obj, cmode = 'outdegree')
  nonsingletons <- which(in_degree != 0 | out_degree != 0)
  return(mean(out_degree[nonsingletons] >= 3))
}
super_spreaders_resample_1 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_before, 
                                     resample_func=resample_by_cluster, 
                                     attribute_func=prop_super_spreaders)
super_spreaders_resample_2 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_after, 
                                     resample_func=resample_by_cluster, 
                                     attribute_func=prop_super_spreaders)
period1_ci_super <- list(percent=c(0.95, 2.50, 97.50, quantile(super_spreaders_resample_1, c(0.025, 0.975))))
period2_ci_super <- list(percent=c(0.95, 2.50, 97.50, quantile(super_spreaders_resample_2, c(0.025, 0.975))))
# CI for out-degree by period
degree_period <- data.table(id = 1:1349, out_degree = out_degree)
degree_period[, period := 2]
degree_period[id %in% (1:1349)[COVID19$origin_time <= TIME_CUT], period := 1]
degree_period <- degree_period[id %in% nonsingletons]
mean(degree_period$out_degree)
out_degree_resample <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net, 
                              resample_func=resample_by_cluster, 
                              attribute_func=function(x){sna::degree(x, cmode='outdegree')})
quantile(out_degree_resample, c(0.025, 0.975))
out_degree_resample_1 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_before, 
                                resample_func=resample_by_cluster, 
                                attribute_func=function(x){sna::degree(x, cmode='outdegree')})
out_degree_resample_2 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_after, 
                                resample_func=resample_by_cluster, 
                                attribute_func=function(x){sna::degree(x, cmode='outdegree')})
period1_ci_degree <- list(percent=c(0.95, 2.50, 97.50, quantile(out_degree_resample_1, c(0.025, 0.975))))
period2_ci_degree <- list(percent=c(0.95, 2.50, 97.50, quantile(out_degree_resample_2, c(0.025, 0.975))))
# CI for diameter by period
feature <- 'diameter'
# overall interval
mean(summary_graph_feature(net, feature = feature))
diameter_resample <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net, 
                            resample_func=resample_by_cluster, 
                            attribute_func=function(x){summary_graph_feature(x, feature = feature)})
quantile(diameter_resample, c(0.025, 0.975))
diameter_resample_1 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_before, 
                              resample_func=resample_by_cluster, 
                              attribute_func=function(x){summary_graph_feature(x, feature = feature)})
diameter_resample_2 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_after, 
                              resample_func=resample_by_cluster, 
                              attribute_func=function(x){summary_graph_feature(x, feature = feature)})
period1_ci_diameter <- list(percent=c(0.95, 2.50, 97.50, quantile(diameter_resample_1, c(0.025, 0.975))))
period2_ci_diameter <- list(percent=c(0.95, 2.50, 97.50, quantile(diameter_resample_2, c(0.025, 0.975))))

# CI for distance by period
feature <- 'distance'
mean(summary_graph_feature(net, feature = feature))
distance_resample <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net, 
                            resample_func=resample_by_cluster, 
                            attribute_func=function(x){summary_graph_feature(x, feature = feature)})
quantile(distance_resample, c(0.025, 0.975))
distance_resample_1 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_before, 
                              resample_func=resample_by_cluster, 
                              attribute_func=function(x){summary_graph_feature(x, feature = feature)})
distance_resample_2 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_after, 
                              resample_func=resample_by_cluster, 
                              attribute_func=function(x){summary_graph_feature(x, feature = feature)})
period1_ci_distance <- list(percent=c(0.95, 2.50, 97.50, quantile(distance_resample_1, c(0.025, 0.975))))
period2_ci_distance <- list(percent=c(0.95, 2.50, 97.50, quantile(distance_resample_2, c(0.025, 0.975))))

# CI for cluster sizes by period
feature <- 'cluster size'
mean(summary_graph_feature(net, feature = feature))
cluster_size_resample <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net, 
                                resample_func=resample_by_cluster, 
                                attribute_func=function(x){summary_graph_feature(x, feature = feature)})
quantile(cluster_size_resample, c(0.025, 0.975))
cluster_size_resample_1 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_before, 
                                  resample_func=resample_by_cluster, 
                                  attribute_func=function(x){summary_graph_feature(x, feature = feature)})
cluster_size_resample_2 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_after, 
                                  resample_func=resample_by_cluster, 
                                  attribute_func=function(x){summary_graph_feature(x, feature = feature)})
period1_ci_size <- list(percent=c(0.95, 2.50, 97.50, quantile(cluster_size_resample_1, c(0.025, 0.975))))
period2_ci_size <- list(percent=c(0.95, 2.50, 97.50, quantile(cluster_size_resample_2, c(0.025, 0.975))))

# CI for betweenness by period
feature <- 'betweenness'
mean(summary_graph_feature(net, feature = feature))
betweenness_resample <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net, 
                               resample_func=resample_by_cluster, 
                               attribute_func=function(x){summary_graph_feature(x, feature = feature)})
quantile(betweenness_resample, c(0.025, 0.975))
betweenness_resample_1 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_before, 
                                 resample_func=resample_by_cluster, 
                                 attribute_func=function(x){summary_graph_feature(x, feature = feature)})
betweenness_resample_2 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_after, 
                                 resample_func=resample_by_cluster, 
                                 attribute_func=function(x){summary_graph_feature(x, feature = feature)})
period1_ci_betw <- list(percent=c(0.95, 2.50, 97.50, quantile(betweenness_resample_1, c(0.025, 0.975))))
period2_ci_betw <- list(percent=c(0.95, 2.50, 97.50, quantile(betweenness_resample_2, c(0.025, 0.975))))


# overall interval for household
get_household_prop <- function(net_obj){
  household_table <- table(get.edge.attribute(net_obj, 'household'))
  return(household_table['#e41a1c'] / sum(household_table))
}
get_household_prop(net)
household_resample <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net, 
                             resample_func=resample_by_cluster, 
                             attribute_func=get_household_prop)
quantile(household_resample, c(0.025, 0.975))

# CI for household proportion by period
household_resample_1 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_before, 
                               resample_func=resample_by_cluster, 
                               attribute_func=get_household_prop)
household_resample_2 <- sapply(2020 + 1:N_TIMES, bootMean2, net_obj=net_after, 
                               resample_func=resample_by_cluster, 
                               attribute_func=get_household_prop)
period1_ci_house <- list(percent=c(0.95, 2.50, 97.50, quantile(household_resample_1, c(0.025, 0.975))))
period2_ci_house <- list(percent=c(0.95, 2.50, 97.50, quantile(household_resample_2, c(0.025, 0.975))))


# test 
N_TIMES <- 100
stat_super_array <- rep(0, N_TIMES)
stat_t_degree_array <- rep(0, N_TIMES)
stat_degree_array <- rep(0, N_TIMES)
stat_t_diameter_array <- rep(0, N_TIMES)
stat_diameter_array <- rep(0, N_TIMES)
stat_t_dist_array <- rep(0, N_TIMES)
stat_dist_array <- rep(0, N_TIMES)
stat_t_size_array <- rep(0, N_TIMES)
stat_size_array <- rep(0, N_TIMES)
stat_t_betw_array <- rep(0, N_TIMES)
stat_betw_array <- rep(0, N_TIMES)
stat_t_house_array <- rep(0, N_TIMES)

# test super spreaders
super_spreaders <- is_singletons %>%
  subset(is_single == FALSE) %>%
  group_by(period) %>%
  summarise(super = sum(out_degree >= 3),
            n_cases = n(),
            prop = super / n_cases) %>%
  as.data.table()
test <- prop.test(super_spreaders$super, super_spreaders$n_cases, alternative = 'greater')
prop_super <- test$estimate
stat_super <- test$statistic

# test degree
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

test <- t.test(degree_period[period == 1, out_degree],
               degree_period[period == 2, out_degree], alternative = 'greater')
stat_t_degree <- test$statistic
mean_outdegree <- test$estimate

test <- chisq.test(matrix(degree_dist$N, nrow = 2))
stat_degree <- test$statistic

# test diameter
feature <- 'diameter'
diam_period <- rbind(data.table(period = 1, diam = summary_graph_feature(net_before, feature = feature)),
                     data.table(period = 2, diam = summary_graph_feature(net_after, feature = feature)))

diam_dist <- diam_period %>%
  group_by(period) %>%
  summarise('1' = sum(diam == 1),
            '2' = sum(diam == 2),
            '3' = sum(diam == 3),
            '4' = sum(diam == 4)) %>%
  melt(id.vars = "period",
       value.name = "N",
       variable.name = "diam") %>%
  as.data.table()
diam_dist[, total := sum(N), by = period][, prop := N / total][, total := NULL]
diam_dist[, period := factor(period, levels = c(1, 2), labels = c("I", "II"))]

test <- t.test(diam_period[period == 1, diam],
               diam_period[period == 2, diam], alternative = 'greater')
mean_diameter <- test$estimate
stat_t_diameter <- test$statistic
test <- chisq.test(matrix(diam_dist$N, nrow = 2))
stat_diameter <- test$statistic

# test distance
feature <- 'distance'
dist_period <- rbind(data.table(period = 1, dist = summary_graph_feature(net_before, feature = feature)),
                     data.table(period = 2, dist = summary_graph_feature(net_after, feature = feature)))

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

test <- t.test(dist_period[period == 1, dist],
               dist_period[period == 2, dist], alternative = 'greater')
mean_dist <- test$estimate
stat_t_dist <- test$statistic
test <- chisq.test(matrix(dist_dist$N, nrow = 2))
stat_dist <- test$statistic

# test size
feature <- 'cluster size'
size_period <- rbind(data.table(period = 1, size = summary_graph_feature(net_before, feature = feature)),
                     data.table(period = 2, size = summary_graph_feature(net_after, feature = feature)))

size_dist <- size_period %>%
  group_by(period) %>%
  summarise('2' = sum(size == 2),
            '3' = sum(size == 3),
            '4' = sum(size == 4),
            '>=5' = sum(size >= 5)) %>%
  melt(id.vars = "period",
       value.name = "N",
       variable.name = "size") %>%
  as.data.table()
size_dist[, total := sum(N), by = period][, prop := N / total][, total := NULL]
size_dist[, period := factor(period, levels = c(1, 2), labels = c("I", "II"))]

test <- t.test(size_period[period == 1, size],
               size_period[period == 2, size], alternative = 'greater')
mean_size <- test$estimate
stat_t_size <- test$statistic
test <- chisq.test(matrix(size_dist$N, nrow = 2))
stat_size <- test$statistic

# test betw
feature <- 'betweenness'
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

test <- t.test(betw_period[period == 1, betw],
               betw_period[period == 2, betw], alternative = 'greater')
mean_betw <- test$estimate
stat_t_betw <- test$statistic
test <- chisq.test(matrix(betw_dist$N, nrow = 2))
stat_betw <- test$statistic

# test household
household_before <- table(network::get.edge.attribute(net_before, 'household'))
household_after <- table(network::get.edge.attribute(net_after, 'household'))
test <- prop.test(x = c(household_before['#e41a1c'], household_after['#e41a1c']),
                  n = c(sum(household_before), sum(household_after)), alternative = "less")
mean_house <- test$estimate
stat_t_house <- test$statistic



for (seed in 1:N_TIMES){
  tmp <- resample_by_cluster_pool(cluster_before, cluster_after, seed=2000 + seed)
  
  # test super spreaders
  super_before <- get_number_super_spreaders(tmp[[1]])
  super_after <- get_number_super_spreaders(tmp[[2]])
  test <- prop.test(c(super_before[1], super_after[1]), c(super_before[2], super_after[2]), alternative = 'greater')
  stat_super_array[seed] <- test$statistic
  # test out-degree
  out_degree_resample_1 <- sna::degree(tmp[[1]], cmode='outdegree')
  out_degree_resample_2 <- sna::degree(tmp[[2]], cmode='outdegree')
  
  test <- t.test(out_degree_resample_1, out_degree_resample_2, alternative = 'greater')
  stat_t_degree_array[seed] <- test$statistic
  
  degree_period_boot <- rbind(data.table(period=1, out_degree=out_degree_resample_1),
                              data.table(period=2, out_degree=out_degree_resample_2))
  degree_dist_boot <- degree_period_boot %>%
    group_by(period) %>%
    summarise('0' = sum(out_degree == 0),
              '1' = sum(out_degree == 1),
              '2' = sum(out_degree == 2),
              '>=3' = sum(out_degree >= 3)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "degree") %>%
    as.data.table()
  con_table <- matrix(degree_dist_boot$N, nrow=2)
  test <- chisq.test(con_table[, colSums(con_table)>0])
  stat_degree_array[seed] <- test$statistic
  # test diameter
  feature <- 'diameter'
  diameter_resample_1 <- summary_graph_feature(tmp[[1]], feature = feature)
  diameter_resample_2 <- summary_graph_feature(tmp[[2]], feature = feature)
  
  test <- t.test(diameter_resample_1, diameter_resample_2, alternative = 'greater')
  stat_t_diameter_array[seed] <- test$statistic
  
  diam_period_boot <- rbind(data.table(period = 1, diam = diameter_resample_1),
                            data.table(period = 2, diam = diameter_resample_2))
  diam_dist_boot <- diam_period_boot %>%
    group_by(period) %>%
    summarise('1' = sum(diam == 1),
              '2' = sum(diam == 2),
              '3' = sum(diam == 3),
              '4' = sum(diam == 4)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "diam") %>%
    as.data.table()
  con_table <- matrix(diam_dist_boot$N, nrow=2)
  test <- chisq.test(con_table[, colSums(con_table)>0])
  stat_diameter_array[seed] <- test$statistic
  # test distance
  feature <- 'distance'
  dist_resample_1 <- summary_graph_feature(tmp[[1]], feature = feature)
  dist_resample_2 <- summary_graph_feature(tmp[[2]], feature = feature)
  
  test <- t.test(dist_resample_1, dist_resample_2, alternative = 'greater')
  stat_t_dist_array[seed] <- test$statistic
  
  dist_period_boot <- rbind(data.table(period = 1, dist = dist_resample_1),
                            data.table(period = 2, dist = dist_resample_2))
  dist_dist_boot <- dist_period_boot %>%
    group_by(period) %>%
    summarise('1' = sum(dist == 1),
              '2' = sum(dist == 2),
              '3' = sum(dist == 3),
              '4' = sum(dist == 4)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "dist") %>%
    as.data.table()
  con_table <- matrix(dist_dist_boot$N, nrow=2)
  test <- chisq.test(con_table[, colSums(con_table)>0])
  stat_dist_array[seed] <- test$statistic
  # test size
  feature <- 'cluster size'
  size_resample_1 <- summary_graph_feature(tmp[[1]], feature = feature)
  size_resample_2 <- summary_graph_feature(tmp[[2]], feature = feature)
  
  test <- t.test(size_resample_1, size_resample_2, alternative = 'greater')
  stat_t_size_array[seed] <- test$statistic
  
  size_period_boot <- rbind(data.table(period = 1, size = size_resample_1),
                            data.table(period = 2, size = size_resample_2))
  size_dist_boot <- size_period_boot %>%
    group_by(period) %>%
    summarise('2' = sum(size == 2),
              '3' = sum(size == 3),
              '4' = sum(size == 4),
              '>=5' = sum(size >= 5)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "size") %>%
    as.data.table()
  con_table <- matrix(size_dist_boot$N, nrow=2)
  test <- chisq.test(con_table[, colSums(con_table)>0])
  stat_size_array[seed] <- test$statistic
  # test betw
  feature <- 'betweenness'
  betw_resample_1 <- summary_graph_feature(tmp[[1]], feature = feature)
  betw_resample_2 <- summary_graph_feature(tmp[[2]], feature = feature)
  
  test <- t.test(betw_resample_1, betw_resample_2, alternative = 'greater')
  stat_t_betw_array[seed] <- test$statistic
  
  betw_period_boot <- rbind(data.table(period = 1, betw = betw_resample_1),
                            data.table(period = 2, betw = betw_resample_2))
  betw_dist_boot <- betw_period_boot %>%
    group_by(period) %>%
    summarise('2' = sum(betw == 2),
              '3' = sum(betw == 3),
              '4' = sum(betw == 4),
              '>=5' = sum(betw >= 5)) %>%
    melt(id.vars = "period",
         value.name = "N",
         variable.name = "betw") %>%
    as.data.table()
  con_table <- matrix(betw_dist_boot$N, nrow=2)
  test <- chisq.test(con_table[, colSums(con_table)>0])
  stat_betw_array[seed] <- test$statistic
  
  # test household
  household_resample_1 <- table(network::get.edge.attribute(tmp[[1]], 'household'))
  household_resample_2 <- table(network::get.edge.attribute(tmp[[2]], 'household'))
  test <- prop.test(x = c(household_resample_1['#e41a1c'], household_resample_2['#e41a1c']),
                    n = c(sum(household_resample_1), sum(household_resample_2)), alternative = "less")
  stat_t_house_array[seed] <- test$statistic
}

#
if (USE_BOOT){
  test <- list(p.value=mean(stat_super_array >= stat_super), estimate = prop_super)
  test_results <- rbind(test_results, test_to_df("prop_super_spreaders", test, period1_ci_super, period2_ci_super))
  test <- list(p.value=mean(stat_t_degree_array >= stat_t_degree), estimate = mean_outdegree)
  test_results <- rbind(test_results, test_to_df("avg_outdegree", test, period1_ci_degree, period2_ci_degree))
  test <- list(p.value=mean(stat_t_diameter_array >= stat_t_diameter), estimate = mean_diameter)
  test_results <- rbind(test_results, test_to_df("avg_diameter", test, period1_ci_diameter, period2_ci_diameter))
  test <- list(p.value=mean(stat_t_dist_array >= stat_t_dist), estimate = mean_dist)
  test_results <- rbind(test_results, test_to_df("avg_distance", test, period1_ci_distance, period2_ci_distance))
  test <- list(p.value=mean(stat_t_size_array >= stat_t_size), estimate = mean_size)
  test_results <- rbind(test_results, test_to_df("avg_cluster_size", test, period1_ci_size, period2_ci_size))
  test <- list(p.value=mean(stat_t_betw_array >= stat_t_betw), estimate = mean_betw)
  test_results <- rbind(test_results, test_to_df("avg_betweenness", test, period1_ci_betw, period2_ci_betw))
  
  test <- list(p.value=mean(stat_degree_array >= stat_degree))
  test_results <- rbind(test_results, test_to_df("outdegree", test, NULL, NULL))
  test <- list(p.value=mean(stat_diameter_array >= stat_diameter))
  test_results <- rbind(test_results, test_to_df("diameter", test, NULL, NULL))
  test <- list(p.value=mean(stat_dist_array >= stat_dist))
  test_results <- rbind(test_results, test_to_df("distance", test, NULL, NULL))
  test <- list(p.value=mean(stat_size_array >= stat_size))
  test_results <- rbind(test_results, test_to_df("cluster_size", test, NULL, NULL))
  test <- list(p.value=mean(stat_betw_array >= stat_betw))
  test_results <- rbind(test_results, test_to_df("betweenness", test, NULL, NULL))
  
  test <- list(p.value=mean(stat_t_house_array >= stat_t_house), estimate = mean_house)
  test_results <- rbind(test_results, test_to_df("prop_household", test, period1_ci_house, period2_ci_house))
}


####################test for removal periods####################
COVID19$remove_period <- as.Date(COVID19$confirm_date) - as.Date(COVID19$onset_date)
remove_period <- COVID19[, .(time_gap = mean(remove_period)), by = onset_date]
setkey(remove_period, 'onset_date')
remove_period$period <- c(rep('I', 16), rep('II', 31))

# write.csv(remove_period, file="data/remove_period.csv")

test <- funtimes::notrend_test(as.numeric(remove_period$time_gap), test = "MK")

test_results <- rbind(test_results, test_to_df("removal_trend", test, NULL, NULL))

remove_period_period <- data.table(id = 1:1349, remove_period = COVID19$remove_period,
                                   is_symp = ifelse(COVID19$case_dignosis2 == 'confirmed_cases', 1, 0))
remove_period_period[, period := 2]
remove_period_period[id %in% (1:1349)[COVID19$onset_date <= TIME_CUT], period := 1]

remove_period_period %>%
  group_by(period) %>%
  summarise(mean_remove = mean(remove_period),
            sd_remove = sd(remove_period) / sqrt(n()),
            n_cases = n())

test <- t.test(data = remove_period_period, remove_period ~ period, alternative = 'greater')

test_results <- rbind(test_results, test_to_df("avg_removal", test, NULL, NULL))

COVID19 %>%
  subset(onset_date >= '2020-02-10') %>%
  summarise(mean_remove=mean(remove_period),
            se_remove=sd(remove_period) / sqrt(n()))

####################multiple adjustment####################
test_results[, p_value := round(p.adjust(p_value, method = "BH"), 4)]

