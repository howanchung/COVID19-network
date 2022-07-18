library(funtimes)

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

# num of clusters
cluster_before <- igraph::decompose(intergraph::asIgraph(net_before))
cluster_before <- cluster_before[sapply(cluster_before, igraph::vcount) > 1]
cluster_after <- igraph::decompose(intergraph::asIgraph(net_after))
cluster_after <- cluster_after[sapply(cluster_after, igraph::vcount) > 1]
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
boot.ci(boot(is_singletons$is_single, bootMean, R = 10000), type = "perc")

singletons <- is_singletons %>%
    group_by(period) %>%
    summarise(is_singletons = sum(is_single),
              n_cases = n()) %>%
    as.data.table()

singletons

# CI for proportion
period1_ci <- boot.ci(boot(is_singletons[period == 1, is_single], bootMean, R = 10000), type = "perc")
period2_ci <- boot.ci(boot(is_singletons[period == 2, is_single], bootMean, R = 10000), type = "perc")

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
period1_ci <- boot.ci(boot(is_singletons[is_single == FALSE & period == 1, out_degree >= 3], bootMean, R = 10000), type = "perc")
period2_ci <- boot.ci(boot(is_singletons[is_single == FALSE & period == 2, out_degree >= 3], bootMean, R = 10000), type = "perc")

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

# overall interval
mean(degree_period$out_degree)
set.seed(2021)
boot.ci(boot(degree_period[, out_degree], bootMean, R = 10000), type = "perc")

# CI for out-degree by period
degree_period %>%
    group_by(period) %>%
    summarise(mean_outdegree = mean(out_degree))
period1_ci <- boot.ci(boot(degree_period[period == 1, out_degree], bootMean, R = 10000), type = "perc")
period2_ci <- boot.ci(boot(degree_period[period == 2, out_degree], bootMean, R = 10000), type = "perc")

test <- t.test(out_degree ~ period, data = degree_period, alternative = 'greater')
# wilcox.test(out_degree ~ period, data = degree_period, exact = FALSE, alternative = 'greater')

test_results <- rbind(test_results, test_to_df("avg_outdegree", test, period1_ci, period2_ci))

# distribution comparison
ks.test(x = degree_period[period == 1, out_degree],
        y = degree_period[period == 2, out_degree],
        alternative = "greater")

degree_period[, c("out_degree", "period")] %>%
    group_by(period) %>%
    table() %>%
    chisq.test()

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

# overall interval
mean(summary_graph_feature(net, feature = 'diameter'))
set.seed(2021)
boot.ci(boot(summary_graph_feature(net, feature = 'diameter'), bootMean, R = 10000), type = "perc")

# CI for diameter by period
period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = 'diameter'), bootMean, R = 10000), type = "perc")
period2_ci <- boot.ci2(boot(summary_graph_feature(net_after, feature = 'diameter'), bootMean, R = 10000), type = "perc")

#
test <- t.test(x = summary_graph_feature(net_before, feature = 'diameter'),
               y = summary_graph_feature(net_after, feature = 'diameter'),
               alternative = 'greater')

test_results <- rbind(test_results, test_to_df("avg_diameter", test, period1_ci, period2_ci))


data.table(diameter = summary_graph_feature(net_before, feature = 'diameter'), period = 1) %>%
    rbind(data.table(diameter = summary_graph_feature(net_after, feature = 'diameter'), period = 2)) %>%
    group_by(period) %>%
    table() %>%
    chisq.test()

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
# overall interval
mean(summary_graph_feature(net, feature = 'distance'))
set.seed(2021)
boot.ci(boot(summary_graph_feature(net, feature = 'distance'), bootMean, R = 10000), type = "perc")

# CI for distance by period
period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = 'distance'), bootMean, R = 10000), type = "perc")
period2_ci <- boot.ci2(boot(summary_graph_feature(net_after, feature = 'distance'), bootMean, R = 10000), type = "perc")

test <- t.test(x = summary_graph_feature(net_before, feature = 'distance'),
               y = summary_graph_feature(net_after, feature = 'distance'),
               alternative = 'greater')

test_results <- rbind(test_results, test_to_df("avg_distance", test, period1_ci, period2_ci))

wilcox.test(x = summary_graph_feature(net_before, feature = 'distance'),
            y = summary_graph_feature(net_after, feature = 'distance'),
            alternative = 'greater')

ks.test(x = summary_graph_feature(net_before, feature = 'distance'),
        y = summary_graph_feature(net_after, feature = 'distance'))


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
# overall interval
mean(summary_graph_feature(net, feature = 'cluster size'))
set.seed(2021)
boot.ci(boot(summary_graph_feature(net, feature = 'cluster size'), bootMean, R = 10000), type = "perc")

# CI for cluster sizes by period
period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = 'cluster size'), bootMean, R = 10000), type = "perc")
period2_ci <- boot.ci(boot(summary_graph_feature(net_after, feature = 'cluster size'), bootMean, R = 10000), type = "perc")

test <- t.test(x = summary_graph_feature(net_before, feature = 'cluster size'),
               y = summary_graph_feature(net_after, feature = 'cluster size'),
               alternative = 'greater')

test_results <- rbind(test_results, test_to_df("avg_cluster_size", test, period1_ci, period2_ci))

wilcox.test(x = summary_graph_feature(net_before, feature = 'cluster size'),
            y = summary_graph_feature(net_after, feature = 'cluster size'),
            alternative = 'greater')

ks.test(x = summary_graph_feature(net_before, feature = 'cluster size'),
        y = summary_graph_feature(net_after, feature = 'cluster size'))


data.table(diameter = summary_graph_feature(net_before, feature = 'cluster size'), period = 1) %>%
    rbind(data.table(diameter = summary_graph_feature(net_after, feature = 'cluster size'), period = 2)) %>%
    group_by(period) %>%
    table() %>%
    # chisq.test()
    fisher.test()

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
mean(summary_graph_feature(net, feature = feature))
set.seed(2021)
boot.ci(boot(summary_graph_feature(net, feature = feature), bootMean, R = 10000), type = "perc")

# CI for cluster sizes by period
period1_ci <- boot.ci(boot(summary_graph_feature(net_before, feature = feature), bootMean, R = 10000), type = "perc")
period2_ci <- boot.ci(boot(summary_graph_feature(net_after, feature = feature), bootMean, R = 10000), type = "perc")

test <- t.test(x = summary_graph_feature(net_before, feature = feature),
               y = summary_graph_feature(net_after, feature = feature))

test_results <- rbind(test_results, test_to_df("avg_betweenness", test, period1_ci, period2_ci))

wilcox.test(x = summary_graph_feature(net_before, feature = feature),
            y = summary_graph_feature(net_after, feature = feature))

ks.test(x = summary_graph_feature(net_before, feature = feature),
        y = summary_graph_feature(net_after, feature = feature))


data.table(betw = summary_graph_feature(net_before, feature = feature), period = 1) %>%
  rbind(data.table(betw = summary_graph_feature(net_after, feature = feature), period = 2)) %>%
  group_by(period) %>%
  table() %>%
  # chisq.test()
  fisher.test()

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
mean(household)
set.seed(2021)
boot.ci(boot(household, bootMean, R = 10000), type = "perc")

# CI for household proportion by period
set.seed(2021)
period1_ci <- boot.ci(boot(household_before$household, bootMean, R = 10000), type = "perc")
set.seed(2021)
period2_ci <- boot.ci(boot(household_after$household, bootMean, R = 10000), type = "perc")

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

