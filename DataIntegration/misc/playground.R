library(pheatmap)
library(dplyr)
library(reshape2)
library(igraph)
library(linkcomm)
library(kernlab)


devtools::source_url("https://raw.githubusercontent.com/scalefreegan/Teaching/master/DataIntegration/scripts/readData.R")
load("~/Documents/Teaching/DataIntegration/lab.rda")

load("/g/steinmetz/brooks/Teaching/data.rda")
load("/g/steinmetz/brooks/Teaching/lab.rda")

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
  plot.igraph(g, edge.width = E(g)$weight/(max(E(g)$weight)/5), vertex.label = NA, vertex.size = 5)
dev.off()

tmp_c = rep(1,length(mc_spcc2))
tmp_c[mc_spcc2==1] = 2
pdf(file="~/Documents/Teaching/DataIntegration/lab_graph_spec.pdf")
  plot.igraph(g, edge.width = E(g)$weight/(max(E(g)$weight)/5), vertex.label = NA, vertex.size = 5, vertex.color=tmp_c, edge.curved=T)
dev.off()

# graph_laplacian
g_lpc = graph.laplacian(g)

pdf(file="/g/steinmetz/brooks/Teaching/lab_lap.pdf")
  pheatmap(tmp2[1:100,1:100])
dev.off()


mc_spcc = specc(as.kernelMatrix(mc), centers=5)
v = as.data.frame(factor(mc_spcc))
rownames(v) = rownames(mc)
pdf("~/Desktop/tmp.pdf")
  plot.igraph(g,vertex.label=NA,vertex.color=v[names(V(g)),])
dev.off()


# small graph - need to redefine n
mc2 = mc[n,n]
#diag(mc2) = 0
#mc_spcc2 = kkmeans(as.kernelMatrix(mc2), centers=3)
mc_spcc2 = specc(as.kernelMatrix(mc2), centers=5)
v = as.data.frame(factor(mc_spcc2))
rownames(v) = rownames(mc2)
pheatmap(mc2,annotation=v)
g_small = graph_from_adjacency_matrix(mc2, mode = "undirected", weighted = TRUE, diag = FALSE)
plot.igraph(g_small,vertex.label=NA)
plot.igraph(g_small,vertex.label=NA,vertex.color=v[names(V(g_small)),])

g_lpc_eigen = eigen(g_lpc)

library(parallel)
options("mc.cores"=10)
mc_spcc_ = mclapply(seq(5,500,1),function(i){
  print(i)
  try({
    tmp = specc(as.kernelMatrix(mc), centers=i)
    v = as.data.frame(factor(tmp))
    rownames(v) = rownames(mc)
    v = v[names(V(g)),]
    return(modularity(g,v,weights=E(g)$weights))
  })
})
modularity = unlist(mc_spcc_)
names(modularity) = as.character(seq(5,500,1))

pdf(file="~/Documents/Teaching/DataIntegration/modularity.pdf")
  plot(names(modularity),modularity,pch=20, ylab="Modularity (Q)", xlab="# Clusters")
  qm = which(modularity==max(modularity))
  points(names(qm),modularity[qm],col="red")
  text(as.numeric(names(qm))+75,modularity[qm],paste0("Max = ",names(qm)," clusters"),col="red")
dev.off()

png(file="~/Documents/Teaching/DataIntegration/modularity.png")
  plot(names(modularity),modularity,pch=20, ylab="Modularity (Q)", xlab="# Clusters")
  qm = which(modularity==max(modularity))
  points(names(qm),modularity[qm],col="red")
  text(as.numeric(names(qm))+75,modularity[qm],paste0("Max = ",names(qm)," clusters"),col="red")
dev.off()

# link communities
mc_edge = as_edgelist(g)
mc_edge = data.frame(node1 = mc_edge[,1], node2 = mc_edge[,2], weight = E(g)$weight)
g_linkcomm = getLinkCommunities(mc_edge, directed = FALSE)

save(ms, mc, g, mc_edge, g_linkcomm, file="~/Documents/Teaching/DataIntegration/lab.rda")


# derive algorithm on spirals data
library(kernlab)
data(spirals)
# affinoity matrix
# from http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html
s <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])

  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A
}

S = do.call(cbind,lapply(seq(1, dim(spirals)[1]),function(i){
  o = sapply(seq(1, dim(spirals)[1]), function(j){
    if (i!=j) {
      return(s(spirals[i,],spirals[j,]))
    } else {
      return(0)
    }
    })
  return(o)
}))

A <- make.affinity(S, 3)
D <- diag(apply(S, 1, sum))
U <- D - A
k   <- 2
evL <- eigen(U, symmetric=TRUE)
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
k = kmeans(Z,2)
plot(spirals,col=k$cluster)



D = matrix(0,nrow=dim(A)[1],ncol=dim(A)[2])
diag(D) = apply(A,1,sum)

L = solve(D^.5) %*% A %*% solve(D^.5)

L_e = eigen(L,symmetric=T)

X = L_e$vectors[,1:2]

Y = t(apply(X,1,function(i){
  i/sum(i^2)^.5
}))

k = kmeans(Y,2)

plot(spirals,col=k$cluster)

x = spirals
distances = as.matrix(dist(x))
W = exp(-distances^2)
G = diag(rowSums(W))
L = G - W
eig = eigen(L)
km2 = kmeans(cbind(eig$vectors[,298],eig$vectors[,299]),centers=2,iter.max=20,nstart=5)
plot(x,asp=1,xlab="",ylab="",col=km2$cluster)
