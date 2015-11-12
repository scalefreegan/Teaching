library(pheatmap)
library(dplyr)
library(reshape2)
library(igraph)
library(linkcomm)

devtools::source_url("https://raw.githubusercontent.com/scalefreegan/Teaching/master/DataIntegration/readData.R")

data = data_reduced

normalizeKernel = function(K) {
  # from Shawe-Taylor & Cristianini's "Kernel Methods for Pattern Analysis", p113
  # original kernel matrix stored in variable K
  # output uses the same variable K
  # D is a diagonal matrix storing the inverse of the norms
  rnames = rownames(K)
  cnames = colnames(K)
  D = diag(1/sqrt(diag(K)))
  K = D %*% K %*% D
  rownames(K) = rnames
  colnames(K) = cnames
  return(K)
}

ms = lapply(colnames(data)[3:7], function(d){
  z = select(data, gene1, gene2, which(colnames(data)==d))
  z[is.na(z[,d]),d] = 0
  g = sort(unique(c(z$gene1, z$gene2)))
  m = matrix(0, nrow = length(g),ncol = length(g), dimnames = list(g,g))
  m[cbind(z$gene1, z$gene2)] = z[,d]
  m[cbind(z$gene2, z$gene1)] = z[,d]

  # calc eigenvalues
  eigen_m = eigen(m)

  # make sure m is a valid kernel by adding
  # make matrix postitive semi-definite
  toadd = ceiling(abs(min(eigen_m$values)))
  diag(m) = diag(m) + toadd
  nm = normalizeKernel(m)
  o = list()
  o$m = m
  o$nm = nm
  return(o)
})

mc = ms[[1]]$nm
for (i in 2:length(ms)) {
  mc = mc + ms[[i]]$nm
}
mc = normalizeKernel(mc)
# some rounding errors, make symmetric
# mc[lower.tri(mc)] = mc[upper.tri(mc)]

g = graph_from_adjacency_matrix(mc, mode = "undirected", weighted = TRUE, diag = FALSE)

pdf(file="~/Documents/Teaching/DataIntegration/lab_graph.pdf")
  plot.igraph(g, vertex.size = 5, vertex.label = NA)
dev.off()

# link communities
mc_edge = as_edgelist(g)
mc_edge = data.frame(node1 = mc_edge[,1], node2 = mc_edge[,2], weight = E(g)$weight)
g_linkcomm = getLinkCommunities(mc_edge, directed = FALSE)

save(ms, mc, g, mc_edge, g_linkcomm, file="~/Documents/Teaching/DataIntegration/lab.rda")
