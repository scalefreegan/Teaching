# from http://www.di.fc.ul.pt/~jpn/r/spectralclustering/spectralclustering.html

# spirals

library(kernlab)
data(spirals)

s <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

make.similarity <- function(my.data, similarity) {
  N <- nrow(my.data)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      if (i!=j) {
        S[i,j] <- similarity(my.data[i,], my.data[j,])
      } else {
        S[i,j] <- 0
      }
    }
  }
  S
}


S <- make.similarity(spirals, s)

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

A <- make.affinity(S, 4)  # use 3 neighboors (includes self)

D <- diag(apply(A, 1, sum))

U <- D - A

k   <- 2
evL <- eigen(U, symmetric=TRUE)
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]

km <- kmeans(Z, centers=k, nstart=5)
plot(spirals, col=km$cluster)
