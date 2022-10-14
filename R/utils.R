# Title     : TODO
# Objective : TODO
# Created by: howanchung
# Created on: 2021/6/13

####################################################################################################
t_test_by_stat <- function(mu, sd, n, alternative = 'two.sided', conf.level=0.95){
  if (equal_var_test(sd, n) <= 1 - conf.level){
    pool_se <- sqrt(sd[1] ** 2 / n[1] + sd[2] ** 2 / n[2])
    df <- pool_se ** 4 / (sd[1] ** 4 / n[1] ** 2 / (n[1] - 1) + sd[2] ** 4 / n[2] ** 2 / (n[2] - 1))
  }else{
    pool_se <- sqrt(((n[1] - 1) * sd[1] ** 2 + (n[2] - 1) * sd[2] ** 2) / (n[1] + n[2] - 2)) * sqrt(sum(1 / n))
    df <- sum(n) - 2
  }
  
  if (pool_se == 0){
    return(0.0)
  }
  t_stat <- (mu[1] - mu[2]) / pool_se
  if (alternative == "two.sided"){
    p_val <- 2 * min(pt(t_stat, df), 1 - pt(t_stat, df))
  }else if (alternative == "greater"){
    p_val <- 1 - pt(t_stat, df)
  }else if (alternative == "less"){
    p_val <- pt(t_stat, df)
  }else{
    stop("Incorrect alternative")
  }
  return(p_val)
}

equal_var_test <- function(sd, n){
  variance <- sd ** 2
  f_stat <-  variance[1] / variance[2]
  tmp <- pf(f_stat, n[1] - 1, n[2] - 1)
  p <- 2 * min(tmp, 1 - tmp)
  return(p)
}

boot_t_test <- function(mu, se, n, alternative = 'two.sided'){
  # mean_x = N(mu[1], se[1])
  # mean_y = N(mu[2], se[2])
  # (mean_x - mean_y) / sqrt(se[1] + se[2]) = N(0, 1) under H0
  pool_se <- sqrt(sum(se ** 2))
  t_stat <- (mu[1] - mu[2]) / pool_se
  df <- sum(n) - 2
  if (alternative == "two.sided"){
    p_val <- 2 * min(pt(t_stat, df), 1 - pt(t_stat, df))
  }else if (alternative == "greater"){
    p_val <- 1 - pt(t_stat, df)
  }else if (alternative == "less"){
    p_val <- pt(t_stat, df)
  }else{
    stop("Incorrect alternative")
  }
  return(list(estimate=mu, p.value=p_val))
}

resample_by_cluster <- function(net_obj, seed=1234){
  set.seed(seed)
  cluster_list <- igraph::decompose(intergraph::asIgraph(net_obj))
  cluster_list <- cluster_list[sapply(cluster_list, igraph::vcount) > 1]
  n_cluster <- length(cluster_list)
  cluster_id <- sort(sample(1:n_cluster, size=n_cluster, replace=TRUE))
  cluster_resample <- cluster_list[cluster_id]
  return(intergraph::asNetwork(igraph::disjoint_union(cluster_resample)))
}

resample_by_cluster_pool <- function(cluster_list1, cluster_list2, seed=1234){
  set.seed(seed)
  cluster_list <- c(cluster_list1, cluster_list2)
  id <- 1:length(cluster_list)
  cluster_resample1 <- cluster_list[sort(sample(id, size=length(cluster_list1), replace=TRUE))]
  cluster_resample2 <- cluster_list[sort(sample(id, size=length(cluster_list2), replace=TRUE))]
  network_resample1 <- intergraph::asNetwork(igraph::disjoint_union(cluster_resample1))
  network_resample2 <- intergraph::asNetwork(igraph::disjoint_union(cluster_resample2))
  return(list(network_resample1, network_resample2))
}

bootMean <- function(data, indices) mean(data[indices])


bootMean2 <- function(seed, net_obj, resample_func, attribute_func){
  net_resample <- resample_func(net_obj, seed=seed)
  return(mean(attribute_func(net_resample)))
}

bootMean3 <- function(seed, net_obj, resample_func, attribute_func){
  net_resample <- resample_func(net_obj, seed=seed)
  return(attribute_func(net_resample))
}

`%fin%` <- function(x, table) {
    fastmatch::fmatch(x, table, nomatch = 0L) > 0L
}

sample_parent <- function(x) {
    potential_parent1 <- x['contact1_id']
    potential_parent2 <- x['contact2_id']
    potential_parent3 <- x['contact3_id']
    if (is.na(potential_parent1)) {
        return(NA)
    }else if (is.na(potential_parent2)) {
        return(potential_parent1)
    }else if (is.na(potential_parent3)) {
        return(sample(c(potential_parent1, potential_parent2), 1))
    }else {
        return(sample(c(potential_parent1, potential_parent2, potential_parent3), 1))
    }
}

p_value_expr <- function(p_value){
    if (p_value < 0.001){
        return (paste("p-value < 0.001"))
    }else{
        return (paste("p-value =", round(p_value, 3)))
    }
}

active_case <- function(date) {
    return(sum((COVID19$onset_date <= date) & (COVID19$confirm_date >= date)))
}

new_infection <- function(date) {
    return(sum((COVID19$onset_date == date)))
}

summary_graph_feature <- function(T_network, feature = 'diameter') {
  singletons <- which(sna::degree(T_network) == 0)
  T_graph <- intergraph::asIgraph(T_network)
  T_graph <- igraph::delete_vertices(T_graph, singletons)
  if (feature == 'diameter') {
    return(sapply(igraph::decompose(T_graph), igraph::diameter))
  }else if (feature == 'cluster size') {
    return(sapply(igraph::decompose(T_graph), igraph::vcount))
  }else if (feature == 'distance') {
    distances <- igraph::distance_table(T_graph, T)$res
    return(rep(1:length(distances), distances))
  }else if (feature == 'betweenness'){
    cluster <- igraph::decompose(T_graph)
    cluster_size <- sapply(cluster, igraph::vcount)
    cluster <- cluster[cluster_size>1]
    between <- lapply(cluster, igraph::betweenness)
    between <- do.call("c", between)
    return(between)
  }else {
    return(FALSE)
  }
}

summary_graph <- function(net, period, scenario) {
    infected_count <- network::network.size(net)
    singletons <- which(sna::degree(net) == 0)
    num_singletons <- length(singletons)
    p_singletons <- round(num_singletons / infected_count, 5)

    T_graph <- intergraph::asIgraph(net)
    T_graph <- igraph::delete_vertices(T_graph, singletons)
    dist_T <- igraph::mean_distance(T_graph, directed=T)
    diameter_T <- sapply(igraph::decompose(T_graph), igraph::diameter)
    mean_diameter <- mean(diameter_T[diameter_T != 0])
    out_degree <- sna::degree(net, cmode = 'outdegree')
    all_degree <- sna::degree(net, cmode = 'freeman')
    mean_outdegree <- mean(out_degree[all_degree != 0])

    cluster_size <- sapply(igraph::decompose(T_graph), igraph::vcount)
    mean_cluster_size <- mean(cluster_size)
    p_household <- as.vector(table(network::get.edge.attribute(net, 'household')))
    p_household <- p_household[1] / sum(p_household)
    median_ <- c(infected_count, num_singletons, p_singletons, mean_outdegree,
                 mean_diameter, mean_cluster_size, dist_T, p_household)
    Measures <- c('Number of infection', 'Number of singletons', 'Percentage of singletons',
                  'Average outdegree', 'Average diameter of clusters', 'Average size of clusters',
                  'Average shortest path length', 'Proportion within household')
    dt <- data.table(day = period, Measures = Measures, min = NA, median = median_, max = NA, scenario = scenario)
    return(dt)
}

creat_network <- function(seed_, all_contact_infect1, id_information1, n = 20000) {
    all_contact_infect_ <- all_contact_infect1[seed == seed_]
    id_information_ <- id_information1[seed == seed_]
    node_age <- rbind(all_contact_infect_[, c('id', 'id_age')],
                      all_contact_infect_[, c('contact', 'contact_age')], use.names = FALSE)
    node_age <- unique(node_age, by = c('id', 'id_age'))
    setkey(node_age, 'id')
    not_id <- setdiff(1:n, node_age$id)
    not_id_age <- sample(1:14, size = length(not_id), replace = TRUE)
    node_age <- rbind(node_age, data.table(id = not_id, id_age = not_id_age))
    setkey(node_age, 'id')

    all_contact_infect_ <- all_contact_infect_[day <= 47]
    all_contact_infect_ <- unique(all_contact_infect_, by = c('id', 'contact'))
    all_contact_infect_[, id_household := plyr::mapvalues(id, from = id_information_$id,
                                                          to = id_information_$household, warn_missing = F)]
    all_contact_infect_[, contact_household := plyr::mapvalues(contact, from = id_information_$id,
                                                               to = id_information_$household, warn_missing = F)]
    all_contact_infect_[, household := ifelse(id_household == contact_household, 1, 0)]

    household_matrix <- sparseMatrix(all_contact_infect_$id,
                                     all_contact_infect_$contact, x = all_contact_infect_$household,
                                     dims = c(n, n))
    adjacency_matrix <- sparseMatrix(all_contact_infect_$id,
                                     all_contact_infect_$contact, x = 1,
                                     dims = c(n, n))
    in_degree <- colSums(adjacency_matrix)
    node_w_more_than_one_parent <- which(in_degree > 1)
    if (length(node_w_more_than_one_parent) > 0) {
        for (j in node_w_more_than_one_parent) {
            parent_j <- which(adjacency_matrix[, j] == 1)
            true_parent_j <- sample(parent_j, size = 1)
            not_parent_j <- setdiff(parent_j, true_parent_j)
            adjacency_matrix[not_parent_j, j] <- 0
            adjacency_matrix[true_parent_j, j] <- 1
        }
    }
    T_graph <- igraph::graph.adjacency(adjacency_matrix, mode = 'directed')
    age_interval <- plyr::mapvalues(node_age$id_age, from = 1:14, to = rep(1:7, each = 2))
    T_network <- intergraph::asNetwork(T_graph)
    network::set.vertex.attribute(T_network, 'age_interval', age_interval)
    network::set.edge.value(T_network, "household",
                            ifelse(household_matrix == 1, "red", "grey75"))
    return(T_network)
}

sum_graph_dt_period_covid <- function(time_end, time_begin, all_contact_infect_seed,
                                      day_information_seed, id_information_seed, n) {
    all_contact_infect1 <- all_contact_infect_seed[(day > time_begin) & (day <= time_end)]
    all_contact_infect1 <- unique(all_contact_infect1, by = c('id', 'contact'))
    if (nrow(all_contact_infect1) > 0) {
        adjacency_matrix <- sparseMatrix(all_contact_infect1$id,
                                         all_contact_infect1$contact, x = 1,
                                         dims = c(n, n))
        in_degree <- colSums(adjacency_matrix)
        node_w_more_than_one_parent <- which(in_degree > 1)
        if (length(node_w_more_than_one_parent) > 0) {
            for (j in node_w_more_than_one_parent) {
                parent_j <- which(adjacency_matrix[, j] == 1)
                true_parent_j <- sample(parent_j, size = 1)
                adjacency_matrix[, j] <- 0
                adjacency_matrix[true_parent_j, j] <- 1
            }
        }
        in_degree <- colSums(adjacency_matrix)
        out_degree <- rowSums(adjacency_matrix)
        all_degree <- in_degree + out_degree
        T_graph <- igraph::graph.adjacency(adjacency_matrix, mode = 'directed')

        infected <- (1:n)[between(id_information_seed$infected_day, time_begin, time_end + 0.1)]
        infected_count <- length(infected)

        singletons <- which(all_degree == 0)
        T_graph <- igraph::delete_vertices(T_graph, singletons)

        singletons <- singletons[singletons %fin% infected]
        num_singletons <- length(singletons)
        p_singletons <- round(num_singletons / infected_count, 5)

        dist_T <- igraph::mean_distance(T_graph, directed=T)
        diameter_T <- sapply(igraph::decompose(T_graph), igraph::diameter)
        mean_diameter <- mean(diameter_T[diameter_T != 0])
        mean_outdegree <- mean(out_degree[all_degree != 0])
        cluster_size <- sapply(igraph::decompose(T_graph), igraph::vcount)
        mean_cluster_size <- mean(cluster_size)
        return(data.table(period = time_end,
                          new_symp = day_information_seed$new_symp[time_end],
                          singletons = num_singletons,
                          p_singletons,
                          mean_outdegree,
                          mean_diameter,
                          mean_cluster_size,
                          dist_T))
    }else {
        return(NULL)
    }
}

sum_graph_dt_period <- function(time_end, time_begin, all_contact_infect_seed,
                                day_information_seed, id_information_seed, n) {
    all_contact_infect1 <- all_contact_infect_seed[(day > time_begin) & (day <= time_end)]
    all_contact_infect1 <- unique(all_contact_infect1, by = c('id', 'contact'))
    if (nrow(all_contact_infect1) > 0) {
        adjacency_matrix <- sparseMatrix(all_contact_infect1$id,
                                         all_contact_infect1$contact, x = 1,
                                         dims = c(n, n))
        in_degree <- colSums(adjacency_matrix)
        node_w_more_than_one_parent <- which(in_degree > 1)
        if (length(node_w_more_than_one_parent) > 0) {
            for (j in node_w_more_than_one_parent) {
                parent_j <- which(adjacency_matrix[, j] == 1)
                true_parent_j <- sample(parent_j, size = 1)
                adjacency_matrix[, j] <- 0
                adjacency_matrix[true_parent_j, j] <- 1
            }
        }
        in_degree <- colSums(adjacency_matrix)
        out_degree <- rowSums(adjacency_matrix)
        all_degree <- in_degree + out_degree
        T_graph <- igraph::graph.adjacency(adjacency_matrix, mode = 'directed')

        infected <- (1:n)[between(id_information_seed$infected_day, time_begin, time_end + 0.1, incbounds = FALSE)]
        infected_count <- length(infected)

        singletons <- which(all_degree == 0)
        T_graph <- igraph::delete_vertices(T_graph, singletons)

        # singletons <- singletons[singletons %fin% infected]
        # num_singletons <- length(singletons)
        # p_singletons <- round(num_singletons / infected_count, 5)

        dist_T <- igraph::mean_distance(T_graph, directed=T)
        diameter_T <- sapply(igraph::decompose(T_graph), igraph::diameter)
        mean_diameter <- mean(diameter_T[diameter_T != 0])
        mean_outdegree <- mean(out_degree[all_degree != 0])
        cluster_size <- sapply(igraph::decompose(T_graph), igraph::vcount)
        mean_cluster_size <- mean(cluster_size)
        return(data.table(period = time_end,
                          infected_count = infected_count / n,
                          new_symp = day_information_seed$new_symp[time_end],
                          # singletons = num_singletons,
                          # p_singletons,
                          mean_outdegree,
                          mean_diameter,
                          mean_cluster_size,
                          dist_T))
    }else {
        return(NULL)
    }
}

summary_graph_period <- function(seed_, all_contact_infect, day_information, id_information, n) {
    all_contact_infect_seed <- all_contact_infect[seed == seed_]
    id_information_seed <- id_information[seed == seed_]
    day_information_seed <- day_information[seed == seed_]
    sum_dt_list <- rbindlist(lapply(1:60, sum_graph_dt_period,
                                    time_begin = 0,
                                    all_contact_infect_seed = all_contact_infect_seed,
                                    day_information_seed = day_information_seed,
                                    id_information_seed = id_information_seed,
                                    n = n))
    sum_dt_list[is.na(new_symp), new_symp := 0]
    setnames(sum_dt_list, c('period', 'Percentage of infection', 'Number of the new-onset', 'Average outdegree',
                            'Average diameter of clusters', 'Average size of clusters',
                            'Average shortest path length'))
    return(sum_dt_list)
}

summary_graph_t_dt <- function(scenario_, simulation_summary) {
    period_ <- unique(simulation_summary$period)
    a <- list()
    for (i in seq_along(period_)) {
        a[[i]] <- data.table(t(apply(simulation_summary[period == period_[i], -1], 2, quantile,
                                     probs = c(0.025, 0.5, 0.975))), keep.rownames = T)
        a[[i]] <- cbind(data.table(day = i), a[[i]])
        setnames(a[[i]], c('day', 'Measures', 'min', 'median', 'max'))
    }
    a <- rbindlist(a)
    a[, scenario := scenario_]
    return(a)
}

household_proportion <- function(seed_, all_contact_infect, id_information) {
    all_contact_infect_seed <- all_contact_infect[seed == seed_]
    id_information_seed <- id_information[seed == seed_]
    all_contact_infect_seed[, id_household
                                  := plyr::mapvalues(id, from = id_information_seed$id,
                                                     to = id_information_seed$household, warn_missing = F)]
    all_contact_infect_seed[, contact_household
                                  := plyr::mapvalues(contact, from = id_information_seed$id,
                                                     to = id_information_seed$household, warn_missing = F)]

# all_contact_infect_seed %>%
#     merge(id_information_seed[, c("id", "household")], by="id", sort=FALSE) %>%
#     merge(id_information_seed[, c("id", "household")], by.x="contact", by.y="id", sort=FALSE) %>%

    household_rate <- sapply(1:60, function(x, trans_pair) {
        trans_pair1 <- trans_pair[day <= x]
        return(mean(trans_pair1$id_household == trans_pair1$contact_household))
    }, all_contact_infect_seed)
    return(household_rate)
}



household_proportion_dt <- function(scenario_, times, all_contact_infect, id_information) {
    household_dt <- do.call('rbind', lapply(1:times, household_proportion, all_contact_infect, id_information))
    household_dt <- data.table(day = 1:60, t(apply(household_dt, 2, quantile,
                                                   probs = c(0.025, 0.5, 0.975), na.rm = T)))
    household_dt[, Measures := 'Proportion within household']
    setnames(household_dt, c('day', 'min', 'median', 'max', 'Measures'))
    setcolorder(household_dt, c('day', 'Measures', 'min', 'median', 'max'))
    household_dt[, scenario := scenario_]
    return(household_dt)
}

symp_curve_new_CI <- function(seed_, day_information) {
    new_symp_curve <- day_information[seed == seed_, c('day', 'new_symp')]
    return(new_symp_curve)
}

symp_curve_new_CI_wrap <- function(day_information) {
    symp_curve_new_CI_list <- lapply(1:100, symp_curve_new_CI, day_information)
    max_date <- max(unlist(lapply(symp_curve_new_CI_list, function(curve) { return(max(curve[, 1])) })))
    new_symp_curve_CI <- data.table(date = 1:max_date, low = rep(0, max_date),
                                    max = rep(0, max_date), mean = rep(0, max_date))
    for (i in 1:max_date) {
        new_confirmed_num <- lapply(symp_curve_new_CI_list,
                                    function(curve) { return(curve[day == i, 2]) })
        new_confirmed_num[sapply(new_confirmed_num, length) == 0] <- 0
        new_confirmed_num <- unlist(new_confirmed_num)
        if (var(new_confirmed_num) == 0) {
            new_symp_curve_CI[i, 2] <- new_symp_curve_CI[i, 3] <- new_symp_curve_CI[i, 4] <- mean(new_confirmed_num)
        }else {
            quantile_curve <- quantile(new_confirmed_num, prob = c(0.025, 0.5, 0.975))
            new_symp_curve_CI[i, 2] <- quantile_curve[1]
            new_symp_curve_CI[i, 3] <- quantile_curve[3]
            new_symp_curve_CI[i, 4] <- quantile_curve[2]
        }
    }
    return(new_symp_curve_CI)
}

symp_curve_count_CI <- function(seed_, day_information) {
    new_symp_curve <- day_information[seed == seed_, c('day', 'count_symp')]
    return(new_symp_curve)
}

symp_curve_count_CI_wrap <- function(day_information) {
    symp_curve_new_CI_list <- lapply(1:100, symp_curve_count_CI, day_information)
    max_date <- max(unlist(lapply(symp_curve_new_CI_list, function(curve) { return(max(curve[, 1])) })))
    new_symp_curve_CI <- data.table(date = 1:max_date, low = rep(0, max_date),
                                    max = rep(0, max_date), mean = rep(0, max_date))
    for (i in 1:max_date) {
        new_confirmed_num <- lapply(symp_curve_new_CI_list,
                                    function(curve) { return(curve[day == i, 2]) })
        new_confirmed_num[sapply(new_confirmed_num, length) == 0] <- 0
        new_confirmed_num <- unlist(new_confirmed_num)
        if (var(new_confirmed_num) == 0) {
            new_symp_curve_CI[i, 2] <- new_symp_curve_CI[i, 3] <- new_symp_curve_CI[i, 4] <- mean(new_confirmed_num)
        }else {
            quantile_curve <- quantile(new_confirmed_num, prob = c(0.025, 0.5, 0.975))
            new_symp_curve_CI[i, 2] <- quantile_curve[1]
            new_symp_curve_CI[i, 3] <- quantile_curve[3]
            new_symp_curve_CI[i, 4] <- quantile_curve[2]
        }
    }
    return(new_symp_curve_CI)
}

estimate_R_user <- function(seed_, day_information, all_contact_infect, id_information, index_case_day) {
    incidence <- day_information[seed == seed_]$new_symp
    all_contact_infect_seed <- all_contact_infect[seed == seed_]
    id_information_seed <- id_information[seed == seed_]
    if (length(incidence) >= 100) {
        incidence <- incidence[1:100]
    }else {
        incidence <- c(incidence, rep(0, 100 - length(incidence)))
    }

    a <- all_contact_infect_seed[day <= 100]
    serial_interval <- plyr::mapvalues(a$contact, from = id_information_seed$id,
                                       to = id_information_seed$infectious_day, warn_missing = F) -
        plyr::mapvalues(a$id, from = id_information_seed$id,
                        to = id_information_seed$infectious_day, warn_missing = F)
    serial_interval <- serial_interval[serial_interval >= 0]
    serial_prob <- c(0, as.numeric(table(serial_interval) / length(serial_interval)))

    incid <- data.table(imported = c(index_case_day$count[1:47], rep(0, 13)),
                        local = pmax(incidence - c(index_case_day$count[1:47], rep(0, 13)), 0))

    R_from_data <- EpiEstim::estimate_R(incid,
                              method = "non_parametric_si",
                              config = EpiEstim::make_config(list(si_distr = serial_prob)))
    return(R_from_data$R$`Mean(R)`)
}

Rt_dt <- function(scenario, times, day_information, all_contact_infect, id_information, wuhan_date) {
    R_dt <- do.call('rbind', lapply(2:times, estimate_R_user,
                                    day_information, all_contact_infect,
                                    id_information, wuhan_date))
    R_dt <- data.table(day = 8:100, t(apply(R_dt, 2, quantile,
                                           probs = c(0.025, 0.5, 0.975), na.rm = T)))
    R_dt[, Measures := 'Effective reproduction numbers']
    R_dt[, situation := scenario]
    setnames(R_dt, c('day', 'min', 'median', 'max', 'Measures', 'scenario'))
    setcolorder(R_dt, c('day', 'Measures', 'min', 'median', 'max', 'scenario'))
    return(R_dt)
}
