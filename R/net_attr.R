centr_Betw <- function(graph, cmode='directed'){
    sna::centralization(net, sna::betweenness, cmode=cmode)
}

Closeness <- function(net){
  graph <- intergraph::asIgraph(net)
  graph_list <- igraph::decompose(graph)
  size <- lapply(graph_list, igraph::vcount)
  graph_list <- graph_list[which(size > 1)]
  close <- Reduce("c", lapply(graph_list, igraph::closeness, normalize=TRUE))
  return(mean(close))
}


library(igraph)
library(intergraph)
library(data.table)
net1 <- graph_from_literal(A-+B:C:D:E, B-+F, C-+G)
net2 <- graph_from_literal(A-+B:C:D, B-+E, C-+F, D-+G:H, H-+J)
net3 <- graph_from_literal(A-+B, B-+C:D, C-+E:F, D-+G:H:J)
net4 <- graph_from_literal(A-+B:C, B-+D:E:F:G:H:J, F-+K:L, L-+M, M-+N, 
                           G-+O:P, O-+Q, Q-+R:S:T)

net_list <- list(net1, net2, net3, net4)

get_net_attr <- function(g, decompose=FALSE){
  out_degree <- mean(degree(g, mode='out'))
  aspl <- mean_distance(g)
  if (decompose){
    cluster <- decompose(g)
    cluster_size <- sapply(cluster, vcount)
    size <- mean(cluster_size)
    cluster <- cluster[cluster_size>1]
    diameter_T <- sapply(cluster, diameter)
    diam <- mean(diameter_T[diameter_T != 0])
    between <- lapply(cluster, igraph::betweenness)
    between <- do.call("c", between)
  }else{
    diam <- diameter(g)
    size <- vcount(g)
    between <- igraph::betweenness(g)
  }
  # between <- between / sum(between)
  between <- mean(between)
  tmp <- list(avg_out_degree=out_degree, 
              avg_shortest_path_length=aspl,
              betweenness=between,
              diameter=diam,
              cluster_size=size)
  return(tmp)
}

get_net_attr(net1)
get_net_attr(net2)
get_net_attr(net3)
get_net_attr(net4)

data.table(do.call('rbind', lapply(net_list, get_net_attr)))


get_net_attr(asIgraph(net), TRUE, betw_cmode='directed')
get_net_attr(asIgraph(net_before), TRUE, betw_cmode='directed')
get_net_attr(asIgraph(net_after), TRUE, betw_cmode='directed')

data.table(do.call('rbind', lapply(list(asIgraph(net), 
                                        asIgraph(net_before), 
                                        asIgraph(net_after)), 
                                   get_net_attr, 
                                   decompose=TRUE)))

mean(harmonic_centrality(asIgraph(net)))
mean(harmonic_centrality(asIgraph(net_before)))
mean(harmonic_centrality(asIgraph(net_after)))

net_cluster <- decompose(asIgraph(net))
net_cluster_vcount <- sapply(net_cluster, vcount)
net_cluster <- net_cluster[net_cluster_vcount>2]
between <- do.call('c', 
                   lapply(net_cluster,
                          function(x){sna::betweenness(asNetwork(x), 
                                                       cmode='directed')}))


