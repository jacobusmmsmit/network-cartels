# install.packages("devtools")
# library("devtools")
# install_github("galanisl/NetHypGeom")
library("NetHypGeom")
library(dplyr)
library(igraph)
net <- ps_model(N = 12, avg.k = 3, gma = 2.3, Temp = 0.15)
# plot_degree_distr(net$network)
# compare_clustering_to_er(net$network, 100, "PS network")
# plot_hyperbolic_net(network = net$network, nodes = net$polar, node.colour = net$polar$theta)
# conn <- get_conn_probs(net = net$network, polar = net$polar, bins = 15)
# plot(conn$dist, conn$prob, pch = 16, xlab = "Hyperbolic distance", ylab = "Connection probability")
# 
# # To embed the network using HyperMap, we set LaBNE+HM's window to 2*pi
hm <- labne_hm(net = net$network, gma = 2.3, Temp = 0.15, k.speedup = 1, w = 2*pi)
# 
# # To embed with LaBNE+HM, we reduce HyperMap's search space from 2*pi 
# # to a small window of 15 degrees around LaBNE's angles
# lh <- labne_hm(net = net$network, gma = 2.3, Temp = 0.15, k.speedup = 10, w = pi/12)
# 
# # Comparison between real and HyperMap-inferred angles and real and LaBNE+HM-inferred angles
# plot(net$polar$theta, hm$polar$theta, pch = 16, 
#      xlab = "Real angles", ylab = "Inferred angles", main = "HyperMap")
# plot(net$polar$theta, lh$polar$theta, pch = 16, 
#      xlab = "Real angles", ylab = "Inferred angles", main = "LaBNE+HM")

g = hm$network
# g = make_full_graph(5, directed = T) %>%
#     add_vertices(5) %>%
#     add_edges(c(5,6, 6,7, 6,8, 8,9, 9,10, 7,10)) %>%
#     add_edges(c(10,9, 9,8, 8,6, 7,6, 6,5))
plot.igraph(g)
write_graph(g, "data/toy_network.txt", "edgelist")
lh <- labne_hm(net = g, gma = 2.3, Temp = 0.15, k.speedup = 10, w = pi/12)
write.csv(lh$polar, file="data/toy_network_polar.csv", row.names = F)

toy_BR_el = read.csv("out/toy_BR_network.txt", sep=",", skip = 1, header = F, col.names=c("src", "dest")) %>% data.matrix
toy_BR = graph_from_edgelist(toy_BR_el)
plot.igraph(toy_BR)

p2p_el <- read.table("data/p2p-Gnutella31.txt")
p2p <- graph.data.frame(p2p_el, directed=T)
hm <- labne_hm(net = p2p, gma = 2.3, Temp = 0.15, k.speedup = 1, w = 2*pi)
