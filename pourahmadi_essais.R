library(randcorr)
library(igraph)
library(gmat)

randcorr_force_zeros <- function(m)
{
  p <- nrow(m)
  
  # Check inputs
  if (length(p)>1 || p<2 || p%%1)
  {
    stop("p must be a scalar integer greater than 1.")
  }
  
  # Step 1 - generate angles theta from PDF (sin(theta))^k, k>=1, 0<theta<pi
  e = matrix(1,p,1)
  theta = matrix(0,p,p)
  for(j in 1 : (p-1)) {
    theta[(j+1):p,j] = randcorr.sample.sink( (p-j)*e[(j+1):p] )
  }
  
  # Sets 0's according to the adjacency matrix
  for (i in 1:nrow(m)) {
    for (j in 1:i) {
      
      # Check if the value is 0
      if (m[i, j] == 0 && i!=j) {
        
        # Impose a pi/2 in the theta matrix in the same position
        theta[i,j] = pi/2
      }
    }
  }

  # Step 2 - construct lower triangular Cholesky factor
  L = matrix(1,p,p)
  for (i in 2 : p) {
    L[i,2:i] = cumprod( sin(theta[i,1:(i-1)]) )
  }
  
  R = cos(theta)
  R[upper.tri(R)] = 0
  
  L = L * R
  
  # Form correlation matrix
  C = L %*% t(L)
  
  return(C)
}

# An example of a graph which does not work
############################################
adj <- matrix(0,7,7)
edges <- rbind(c(1,2),c(1,3),c(2,3),c(2,4),c(2,5),c(3,6),c(3,7),c(4,5),c(6,7))
adj[edges] <- 1
adj <- adj+t(adj)
g.adj <- graph_from_adjacency_matrix(adj)
is.chordal(g.adj)$chordal
plot(as.undirected(g.adj))
plot(as.undirected(ug_to_dag(g.adj)),vertex.label.color = "black",
     vertex.label.cex=2.5,vertex.size=20, vertex.color="lightblue", edge.color="black") # ug_to_diag does not change anything
diag(adj) <- 1

# generation of the adjacency matrix with Pourhamadi, 
# forcing the zeroes in the cholesky matrix
adj_gen <- randcorr_force_zeros(adj)

#resulting adjacency
print(1*(zapsmall(adj_gen, digits = 15)!=0))
# wanted adjacency
adj
# number of additional edges
sum(1*(zapsmall(adj_gen, digits = 15)!=0)!=adj)/2



# An example of a graph which works
############################################
adj <- matrix(0,7,7)
edges <- rbind(c(1,2),c(1,3),c(2,3),c(3,4),c(3,5),c(4,5),c(5,6),c(5,7),c(6,7))
adj[edges] <- 1
adj <- adj+t(adj)
g.adj <- graph_from_adjacency_matrix(adj)
is.chordal(g.adj)$chordal
plot(as.undirected(g.adj))
plot(as.undirected(ug_to_dag(g.adj)),vertex.label.color = "black",
     vertex.label.cex=2.5,vertex.size=20, vertex.color="lightblue", edge.color="black")
     # ug_to_diag does not change anything
diag(adj) <- 1

# generation of the adjacency matrix with Pourhamadi, 
# forcing the zeroes in the cholesky matrix
adj_gen <- randcorr_force_zeros(adj)

#resulting adjacency
print(1*(zapsmall(adj_gen, digits = 15)!=0))
# wanted adjacency
adj
# number of additional edges
sum(1*(zapsmall(adj_gen, digits = 15)!=0)!=adj)/2

###################################################################################
###################################################################################

## plot of the elliptope
##########################

adj <- diag(3)
edges <- rbind(c(1,3),c(2,3))
edges <- rbind(c(1,2),c(2,3))
adj[edges] <- 1
adj <- adj+t(adj)
diag(adj) <- 0
plot(as.undirected(graph_from_adjacency_matrix(adj)),vertex.label.color = "black",
     vertex.label.cex=1.5,vertex.size=30, vertex.color="lightblue", edge.color="black")

diag(adj) <- 1

## Plot of an elliptope
matrices <- replicate(5000, randcorr_force_zeros(adj), simplify = FALSE)

elliptope_vectors <- t(sapply(matrices, function(mat) {
  upper_triangular_part <- mat[upper.tri(mat, diag = FALSE)]
  return(upper_triangular_part)
}))
# Convert the correlation matrices to vectors in the elliptope
#elliptope_vectors <- t(sapply(matrices, get_elliptope_vector))

#layout=layout.circle,

# Plot the elliptope using plotly
plot_ly(x = elliptope_vectors[, 1], y = elliptope_vectors[, 2], z = elliptope_vectors[, 3], 
        type = "scatter3d", mode = "markers", marker = list(size = 2)) %>% layout(title = 
              "Elliptope of 5000 Random 3x3 Correlation Matrices associated to a graph")

###################################################################################
###################################################################################

adj <- matrix(0,10,10)
edges <- rbind(c(1,3),c(2,3),c(3,4),c(4,5),c(4,6),c(5,7),c(6,7),c(7,8),c(8,9),c(8,10))
#edges <- rbind(1:9,2:10)
adj[edges] <- 1
adj <- adj+t(adj)
diag(adj) <- 1
adj


adj <- matrix(0,7,7)
edges <- rbind(c(1,2),c(1,3),c(2,3),c(2,4),c(2,5),c(3,6),c(3,7),c(4,5),c(6,7))
#edges <- rbind(c(1,2),c(1,3),c(2,3),c(3,4),c(3,5),c(4,5),c(5,6),c(5,7),c(6,7))
adj[edges] <- 1
adj <- adj+t(adj)
g.adj <- graph_from_adjacency_matrix(adj)
is.chordal(g.adj)
plot(as.undirected(g.adj))
plot(as.undirected(ug_to_dag(g.adj)))
diag(adj) <- 1
adj

adj <- matrix(0,3,3)
edges <- rbind(c(1,2),c(1,3))
edges <- rbind(c(1,3),c(2,3))
adj[edges] <- 1
adj <- adj+t(adj)
diag(adj) <- 1
adj



library(igraph)
library(gmat)
n <- 9
p <- 0.5
# generate graph using Erdos-Renyi method
g <- erdos.renyi.game(n, p, type = "gnp", directed = FALSE)
# change in chordal
chg <- ug_to_dag(g)
plot(chg,edge.arrow.size = .25)
chg <- as.undirected(chg)
plot(chg)
is_chordal(chg)
adj <- as.matrix(get.adjacency(chg))
diag(adj) <- 1

adj_gen <- randcorr_force_zeros(nrow(adj), adj)
print(1*(zapsmall(adj_gen, digits = 15)!=0))
adj
sum(1*(zapsmall(adj_gen, digits = 15)!=0)!=adj)/2

perm <- hclust(as.dist(1-adj))$order
adj_perm <- adj[perm,perm]     
adj_perm_gen <- randcorr_force_zeros(nrow(adj_perm), adj_perm)
print(1*(zapsmall(adj_perm_gen, digits = 15)!=0))
adj_perm
sum(1*(zapsmall(adj_perm_gen, digits = 15)!=0)!=adj_perm)/2
