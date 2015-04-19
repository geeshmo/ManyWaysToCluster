require(animation)
require(MASS)
require(MCMCpack)
require(RColorBrewer)
require(mnormt)
require(ggplot2)
require(gridExtra)


#### Generate data #### 

plot.generated.data <- function(generated.data) {
  
  num.clusters <- generated.data$input.params$num.clusters
  data.points <- generated.data$generated.data$data.points
  z_n.data <-  generated.data$generated.data$points.to.clusters
  cluster.centers <- generated.data$generated.data$cluster.centers
  
  data.points.df <- as.data.frame(data.points)
  colnames(data.points.df) <- c("x.col","y.col")
  z_n.fct <- factor(z_n.data, level = seq(1,max(z_n.data)))
  colorCount <- 22
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  pl1 <- ggplot(cbind(data.points.df,z_n.fct), aes(x=x.col, y=y.col, color=z_n.fct))  + geom_point(shape=1) +
   scale_colour_manual(values = getPalette(colorCount)) + ggtitle("Data points") + theme(aspect.ratio=1)  #scale_colour_brewer(palette="Set1")
  pl1.xrange <- ggplot_build(pl1)$panel$ranges[[1]]$x.range
  pl1.yrange <- ggplot_build(pl1)$panel$ranges[[1]]$y.range
  cls <- as.factor(c(1:num.clusters))
  
  cluster.centers.df <- as.data.frame(cluster.centers)
  colnames(cluster.centers.df) <- c("x.col","y.col")                              
  pl2 <-ggplot(cbind(cluster.centers.df, cls), aes(x=x.col, y=y.col, color=cls)) + geom_point(shape=1) + 
    coord_cartesian(xlim = pl1.xrange, ylim = pl1.yrange) + scale_colour_manual(values = getPalette(colorCount)) + ggtitle("Cluster centers") + #scale_colour_brewer(palette="Set1") 
    theme(aspect.ratio=1)
  
  pl3 <- ggplot(data.points.df, aes(x=x.col, y=y.col))  + geom_point(shape=1) + ggtitle("How the algs see it") + theme(aspect.ratio=1) + coord_cartesian(xlim = pl1.xrange, ylim = pl1.yrange)
  
  #grid.arrange(pl1, pl2, pl3, ncol=2, widths = unit(0.5, "npc") , heights=unit(0.5, "npc") )
  grid.arrange(pl1, pl2, pl3, nrow=1 )
}


generate.clusters <- function(num.clusters, num.data.points, m_0, V_0, alpha, S_0.inv, nu_0, is.dims.uncorellated) {
  
  cluster.centers <-  mvrnorm(num.clusters, m_0, V_0)
  
  ### Probability that a point belongs to a particular cluster
  pi <- rdirichlet(1, alpha)[1,]
  
  ### Find out the number of points in each cluster accroding to pi's
  ### make sure there are no zero points for any of the clusters
  repeat
  {
    num.points.in.clusters <- rmultinom(1, num.data.points, pi)
    num.points.in.clusters <- num.points.in.clusters[,1]
    if (min(num.points.in.clusters) > 0) {
      break
    }
  }
  
  ### Generare corr matrix for the clusters
  R.inv <- rWishart(num.clusters, nu_0, S_0.inv)
    
  R <- list()
  for (i in 1:num.clusters) {
    R[[i]] <- solve(R.inv[,,i])
    if (is.dims.uncorellated) {
      param <- 3*rbeta(1,0.3,0.3)
      R[[i]][1,1] <- param #runif(1, 0, 1)
      R[[i]][2,2] <- 3-param #runif(1, 0, 1)
      R[[i]][1,2] <- 0
      R[[i]][2,1] <- 0
    }
  }  
  
  data.points <-c()
  for (i in 1:num.clusters) {
    points.for.cluster <-  mvrnorm(num.points.in.clusters[i], cluster.centers[i,], R[[i]])
    data.points <- rbind(data.points, points.for.cluster)
  }
  rownames(points.for.cluster) <- NULL
  
  z_n.data <- c()
  for (i in 1:length(num.points.in.clusters)) {
    z_n.data <- c(z_n.data, rep(i, num.points.in.clusters[i]))
  }
  
  # permute points so that algorithms do not take advantage of it
  permutation <- sample(1:num.data.points)
  data.points <- data.points[permutation,]
  z_n.data <- z_n.data[permutation]
  
  #return results
  res <- list()
  
  res$generated.params <- list(pi = pi,  Sigma = R, point.permutation = permutation)
  res$generated.data <- list(cluster.centers = cluster.centers, num.points.in.clusters = num.points.in.clusters, data.points = data.points, points.to.clusters = z_n.data)
  res
}

# Generates data set given the number of clusters and the total number number of datapoints. 'boundary' parameter controls correlation matrices so that clusters are
# easier to fit in plots
generate.data <- function(num.clusters, num.data.points, boundary, is.dims.uncorellated = FALSE) {
  
  # Priors for cluster center point locations
  m_0 <- c(0,0)
  V_0 <- matrix(c(10*boundary,0,0,10*boundary),2,2)
  
  # Priors for cluster probabilities
  alpha <- seq(100, 100 + (num.clusters-1)*50, 50)
  
  # Prior for rwish (for cluster configurations)
  S_0.inv <- toeplitz((2:1)/15)
  
  nu_0 <- 3
  
  input.params <- list(num.clusters = num.clusters, num.data.points = num.data.points, m_0 = m_0, V_0 = V_0, alpha = alpha, S_0.inv = S_0.inv, nu_0 = nu_0)
  
  res <- generate.clusters(num.clusters, num.data.points, m_0, V_0, alpha, S_0.inv, nu_0, is.dims.uncorellated)
  res$input.params <- input.params
  res
}

### Sample usage
# num.clusters <- 5
# num.data.points <- 1000
# boundary <-1
# res <- generate.data(num.clusters, num.data.points, boundary, TRUE)
# plot.generated.data(res)

#saving/reading the data set to/from a file (not used)
save.data.set <- function(folder, file.name, data.set) {
  setwd(folder)
  saveRDS(data.set, file.name) 
}

read.data.set <- function(folder, file.name) {
  setwd(folder)
  res <- readRDS(file.name)
}

#fldr <- "<folder-to-use>"
#save.data.set(fldr, "tree.Rd", res)
#res <- read.data.set(fldr, "for.lac")


# Plots the graphs for either an iteration or the final result (if iter param is NULL). The graphs generated by this function are used as frames in animation.
# x.df is the dataframe of data points, z_n is point to cluster allocation, z_n.data is the true allocation 
plot.cur.results <- function(x.df, z_n, z_n.data, current.cluster.centers.df, actual.cluster.centers.df, title, iter, custom.cluster.colors = NULL) {
  # Plotting
  cur.clusters <- factor(z_n)

  pl1 <- ggplot(cbind(x.df, cur.clusters), aes(x=x.col, y=y.col, color=cur.clusters))  + geom_point(shape=1) + ggtitle("Current clusters")
  actual.clusters <- factor(z_n.data, level = seq(1,max(z_n.data)))
  pl2 <- ggplot(cbind(x.df, actual.clusters), aes(x=x.col, y=y.col, color=actual.clusters))  + geom_point(shape=1) + 
    scale_colour_brewer(palette="Set1") + ggtitle("Actual clusters")
  
  pl1.xrange <- ggplot_build(pl1)$panel$ranges[[1]]$x.range
  pl1.yrange <- ggplot_build(pl1)$panel$ranges[[1]]$y.range
  
  clst.ctr.enm <- as.factor(sort(unique(z_n)))
  act.clst.ctr.enm <- as.factor(sort(unique(z_n.data)))
  
  num.points.in.clusters <- table(z_n) 
  actual.num.points.in.clusters <- table(z_n.data)
  
  pl3 <-ggplot(cbind(current.cluster.centers.df, clst.ctr.enm, num.points.in.clusters), aes(x=x.col, y=y.col, color=clst.ctr.enm)) + 
    geom_point(shape=1, size=pmax(1,round(num.points.in.clusters/100))) + coord_cartesian(xlim = pl1.xrange, ylim = pl1.yrange) + ggtitle("Current cluster centers")
        
  pl4 <- ggplot(cbind(actual.cluster.centers.df, act.clst.ctr.enm, actual.num.points.in.clusters), aes(x=x.col, y=y.col, color=act.clst.ctr.enm)) +
    geom_point(shape=1, size=pmax(1,round(actual.num.points.in.clusters/100))) + coord_cartesian(xlim = pl1.xrange, ylim = pl1.yrange) + scale_colour_brewer(palette="Set1") +
    ggtitle("Actual cluster centers")
  
  if (is.null(custom.cluster.colors)) {
    pl1 <- pl1 + scale_colour_brewer(palette="Set1")
    pl3 <- pl3 + scale_colour_brewer(palette="Set1")
  } else {
    pl1 <- pl1 + scale_colour_manual(values = custom.cluster.colors)
    pl3 <- pl3 + scale_colour_manual(values = custom.cluster.colors)
  }
  g.pl1 <- ggplotGrob(pl1)
  g.pl2 <- ggplotGrob(pl2)
  g.pl3 <- ggplotGrob(pl3)
  g.pl4 <- ggplotGrob(pl4)
  
  #http://stackoverflow.com/questions/19827583/assign-grid-arrange-to-object
  if (!is.null(iter)){
    title <- paste(title, "iter:", iter)
  } 
  arrangeGrob(g.pl2, g.pl1, g.pl4, g.pl3, ncol=2, main=title)
}

# Produces the animation, the function is called from server.R
generate.animation <- function(alg.run.res) {
  cur.dir <- getwd()
  sub.dir <- "www"
  new.dir <- file.path(cur.dir, sub.dir)
  dir.create(new.dir, showWarnings = FALSE)
  setwd(new.dir)
  htmlfile <- "run_animation.html"
  # removes data from the previous runs
  unlink("images/*", recursive=TRUE)
  unlink("css/*", recursive=TRUE)
  unlink("js/*", recursive=TRUE)
  file.remove(htmlfile)
  saveHTML({
    ani.options(interval = 1, nmax = 50)
    par(mar = c(4, 4, .1, 0.1), mgp = c(2, 0.7, 0))
    lapply(alg.run.res$plot.list,  function(x) {grid.arrange(x, ncol=1)})
  }, img.name = "alg_plot", ani.height = 500, ani.width = 600, autoplay=FALSE, autobrowse=FALSE, verbose=FALSE, loop=FALSE, htmlfile=htmlfile,
  title = "Clustering",
  description = c("Clustering animation"))
  setwd(cur.dir)
  htmlfile
}

#extract.results <- function(alg.run.res) {
#  alg.run.res$z_n
#}

#### Block Gibbs ####

# Eq reference are for Machine Learning: A Probabilistic Perspective 
# by Kevin P. Murphy

# Compute x.bar_k
# Eq. 24.16
get.x.bar_k <- function(z_n, x, N_k) {
  res <- list()
  for (i in as.integer(names(N_k))) {
    res[[i]] <- colSums(as.numeric(z_n==i)*x)*1.0/N_k[as.character(i)]
  }
  res
}


# Eq. 24.10
sample.z_n <- function(x, z_n, num.points, mu_k, Sigma_k, pi_k, num.clusters) {
  z.res <- c()
  
  #permute points (actually produces worse results)
  #points.permuted <- sample.int(num.points,num.points)
  #for (pnt in 1:num.points) {
  #  i <- points.permuted[pnt]  
  N_k <- table(z_n)
  for (i in 1:num.points) {
    # find the corresponding x
    cur.x <- x[i,]
    probs <- c()
    # create a distribution to sample from
    for (j in 1:num.clusters) {
      points.in.cluster <- N_k[as.character(j)]
      if (!is.na(points.in.cluster)) {
        probs <- c(probs, pi_k[[j]]*dmnorm(cur.x, mu_k[[j]], (Sigma_k[[j]] + t(Sigma_k[[j]]))/2))  
      } else {
        probs <- c(probs, 0)
      }
    }
    sml <- rmultinom(1, size = 1, prob = probs)
    sml <- sml[,1]
    z.res <- c(z.res, which.max(sml))
  }
  z.res
}

# Eq. 24.11
sample.pi_k <- function(alpha_k, z_n) {
  N_k <- table(z_n)
  existing.clusters <- as.integer(names(N_k))
  params <- alpha_k[existing.clusters] + table(z_n)
  pi <- rdirichlet(1, params)[1,]
  res <- list()
  for (i in 1:length(existing.clusters)) {
     res[[existing.clusters[i]]] <- pi[i]
  }
  res
}


####Sampling cluster means 

# get V_k for cluster i
# Eq. 24.13
get.V_k.inv <- function(Sigma_k.inv, m_0, V_0.inv, i, N_k) { 
  V_k.inv <- V_0.inv + N_k[as.character(i)]*Sigma_k.inv
  V_k.inv
}

# Eq. 24.14
get.m_k <- function(x, Sigma_k.inv, z_n, m_0, V_0.inv, i, N_k, x.bar_k, V_k.inv) {
  m.k <- solve(V_k.inv)%*%(N_k[as.character(i)]*Sigma_k.inv%*%as.matrix(x.bar_k[[i]]) + V_0.inv%*%as.matrix(m_0))
  m.k
}

test.get.m_k <- function(res) {
  #construct parameters from the real data
  x <- res$generated.data$data.points
  Sigma_k <- res$generated.params$Sigma
  m_0 <- res$input.params$m_0
  V_0 <- res$input.params$V_0
  p.counts <- res$generated.data$num.points.in.clusters
  z_n <- res$generated.data$points.to.clusters
  num.clusters <- res$input.params$num.clusters
  N_k <- table(z_n)
  x.bar_k <- get.x.bar_k(z_n, x, N_k)
  V_0.inv <- solve(V_0)
  ms <- list()
  for (i in 1:length(p.counts)) {
    Sigma_k.inv <- solve(Sigma_k[[i]])
    V_k.inv <- get.V_k.inv(Sigma_k.inv, m_0, V_0.inv, i, N_k)
    ms[[i]] <- get.m_k(x, Sigma_k.inv, z_n, m_0, V_0.inv, i, N_k, x.bar_k, V_k.inv)
  }
  ms
}


#out <- test.get.m_k(res)
#should be close to these
#print(res$generated.data$cluster.centers)

# Eq. 24.12
sample.mu_k <- function(x, Sigma_k, z_n, m_0, V_0) {    
  res.mu <- c()
  N_k <- table(z_n)
  x.bar_k <- get.x.bar_k(z_n, x, N_k) 
  V_0.inv <- solve(V_0)
  # go thru all existing clusters
  for (i in as.integer(names(N_k))) {
    Sigma_k.inv <- solve(Sigma_k[[i]])
    V_k.inv <- get.V_k.inv(Sigma_k.inv, m_0, V_0.inv, i, N_k)
    m_k <- get.m_k(x, Sigma_k.inv, z_n, m_0, V_0.inv, i, N_k, x.bar_k, V_k.inv) 
    res.mu[[i]] <- mvrnorm(1, m_k, solve(V_k.inv))
  }
  res.mu
}

########### End sampling cluster means 

########### Sampling cluster variances 

# Eq. 24.18
get.S_k <- function(S_0, z_n, mu_k, x) {
  res.S_k <-list()
  num.all.points <- dim(x)[1]
  # go thru all existing clusters
  for (i in as.integer(names(table(z_n)))) {
    m <- S_0 
    for (j in 1:num.all.points) {
      m = m  + as.numeric(z_n==i)[j]*t(matrix(x[j,] - mu_k[[i]], nrow=1))%*%(matrix(x[j,] - mu_k[[i]], nrow=1))
    }
    res.S_k[[i]] = m
  }
  res.S_k
}

test.get.S_k <- function(res, num.clusters) {
  #construct parameters from the real data
  x <- res$generated.data$data.points
  p.counts <- res$generated.data$num.points.in.clusters
  z_n <- res$generated.data$points.to.clusters
  mu_k.init <- res$generated.data$cluster.centers
  mu_k <- list()
  for (i in 1:num.clusters) {
    mu_k[[i]] <- mu_k.init[i,]
  }
  S_0 <- solve(res$input.params$S_0.inv)
  ms <- get.S_k(S_0, z_n, mu_k, x)
  for (i in 1:length(p.counts)) {
    ms[[i]] = cov2cor(ms[[i]])
  }
  ms
}

### The matrices
### here
#test.get.S_k(res, num.clusters)
# and here
#for (i in 1:num.clusters) {
#  print(cov2cor(res$generated.params$Sigma[[i]]))
#}
# should be very similar

# Eq. 24.17
sample.Sigma_k <- function(S_0, z_n, mu_k, x, nu_0) {
  S_k <- get.S_k(S_0, z_n, mu_k, x)
  N_k <- table(z_n)
  sampled.Sigma_k <- list()
  for (i in as.integer(names(N_k))) {
    Sigma_k.inv <- rWishart(1, nu_0 + N_k[as.character(i)], solve(S_k[[i]]))
    sampled.Sigma_k[[i]] <- solve(Sigma_k.inv[,,1])
  }
  sampled.Sigma_k
}

########### End sampling cluster variances

run.block.gibbs <- function(res, num.clusters, num.iterations, progress.fun = NULL) {
  
  x <- res$generated.data$data.points
  num.points <- res$input.params$num.data.points
  # for plotting
  plot.list <-list()
  x.df <- as.data.frame(x)
  colnames(x.df) <- c("x.col","y.col")
  
  m_0 <- res$input.params$m_0
  nu_0 <- res$input.params$nu_0
  V_0 <-  res$input.params$V_0
  S_0 <- solve(res$input.params$S_0.inv) 
  #num.clusters <- res$input.params$num.clusters
  
  # For plotting
  z_n.data <- res$generated.data$points.to.clusters
  actual.cluster.centers.df <- as.data.frame(res$generated.data$cluster.centers)
  colnames(actual.cluster.centers.df) <- c("x.col","y.col")
  
  ## Generate alpha
  #this alpha has no pripr for z_n
  #alpha_k <- rep(1.0/k, k)
  alpha_k <- res$input.params$alpha
  
  ## Generate a random assignment for z
  p.out <- rep (1.0/num.clusters, num.clusters)
  z_n <- sample(1:num.clusters, num.points, prob = p.out, replace = TRUE) 
  
  mu_k.init <-  mvrnorm(num.clusters, m_0, V_0)
  
  mu_k <- list()
  for (i in 1:num.clusters) {
    mu_k[[i]] <- mu_k.init[i,]
  }
  
  cat("starting Gibbs sampling\n")
  for (iter in 1:num.iterations) {
    
    # graph the iter
    mu_k.iter <- c()
    for (i in 1:num.clusters) {
      mu_k.iter <- rbind(mu_k.iter, mu_k[[i]])
    }
    
    cat("iter ", iter, "\n")
    
    # Plot for the current iteration
    current.cluster.centers.df <- as.data.frame(mu_k.iter)
    colnames(current.cluster.centers.df) <- c("x.col","y.col")                              
    cur.plot <- plot.cur.results(x.df, z_n, z_n.data, current.cluster.centers.df, actual.cluster.centers.df, "Block Gibbs", iter)
    plot.list[[length(plot.list)+1]] <- cur.plot
    if (iter %% 5 == 0 || iter == 1) {
      grid.arrange(cur.plot, ncol=1)
      print(table(z_n))
    }
    
    Sigma_k <- sample.Sigma_k(S_0, z_n, mu_k, x, nu_0)
    pi_k <- sample.pi_k(alpha_k, z_n)
    z_n <- sample.z_n(x, z_n, num.points, mu_k, Sigma_k, pi_k, num.clusters)
    mu_k <- sample.mu_k(x, Sigma_k, z_n, m_0, V_0)
    
    if (!is.null(progress.fun)) {
      progress.fun(1/num.iterations, detail = paste("Doing iter", iter))
    }
    
  }
  
  # Plot the final result
  current.cluster.centers.df <- as.data.frame(mu_k.iter)
  colnames(current.cluster.centers.df) <- c("x.col","y.col")                              
  cur.plot <- plot.cur.results(x.df, z_n, z_n.data, current.cluster.centers.df, actual.cluster.centers.df, "Block Gibbs - final", NULL)
  plot.list[[length(plot.list)+1]] <- cur.plot
  
  list(z_n = z_n, plot.list=plot.list)
}

#run.block.gibbs(res, res$input.params$num.clusters, 80)

#### Collapsed Gibbs #### 

# We are using "Gibbs sampling for fitting finite and infinite Gaussian mixture models" by Herman Kamper, all eq refs are referring to that paper
# It can be found at http://www.kamperh.com/notes/kamper_bayesgmm13.pdf

# calculate for i^th cluster 
get.m_k.S_k_i <- function(S_0, m_0, kappa_0, nu_0, z_n, x, i, precalculated.point.sum, precalculated.S, N_i) {
  num.all.points <- dim(x)[1]
  
  x_bar_i <- precalculated.point.sum[[i]]
  
  x_bar_i <- x_bar_i*1.0/N_i
  
  # equation (7)
  m_k_i <- (kappa_0*m_0+N_i*x_bar_i)/(kappa_0+N_i) 
  
  # equation (8)
  S<- matrix(c(0,0,0,0),nrow=2)
  S <- S + precalculated.S[[i]]
  
  S_k_i <- S_0 + S 
  S_k_i <- S_k_i + kappa_0*t(matrix(m_0, nrow=1))%*%(matrix(m_0, nrow=1))
  S_k_i <- S_k_i - (kappa_0 + N_i)*t(matrix(m_k_i, nrow=1))%*%(matrix(m_k_i, nrow=1))

  res <- list()
  res[[1]] <- m_k_i
  res[[2]] <- S_k_i
  res
}

test.get.m_k.S_k_i <- function(res) {
  #construct parameters from the real data
  x <- res$generated.data$data.points
  #Sigma_k <- res$Sigma
  m_0 <- res$input.params$m_0
  S_0 <- solve(res$input.params$S_0.inv)
  nu0 <- res$input.params$nu_0
  z_n <- res$generated.data$points.to.clusters
  precalculated.S <- precalculate.S.for.all.clusters(x, z_n)
  precalculated.point.sum <- precalculate.point.sum.for.all.clusters(x, z_n)
  
  #FIX ME - robust way to find value for kappa_0
  kappa_0 <-0.05
  ms <- list()
  N_k.table <- table(z_n)
  for (i in 1:num.clusters) {
    N_i <- N_k.table[i]
    #ll <- get.m_k.S_k_i(S_0, m_0, kappa_0, nu_0, z_n, x, i)
    ll <- get.m_k.S_k_i(S_0, m_0, kappa_0, nu_0, z_n, x, i, precalculated.point.sum, precalculated.S, N_i)
    print(ll[[1]])
    print(cov2cor(ll[[2]]))
    print("============")
    ms[[i]] <- ll
  }
  ms
}

### compare printed output
#discard <- test.get.m_k.S_k_i(res);
##### with this
#for (i in 1:num.clusters) {
  #print(cov2cor(res$generated.params$Sigma[[i]]))
#}
### and this
# res$generated.data$cluster.centers


####

# equation (26)
# k is the cluster
get.P_z_k <- function(alpha, z_n.minus.i, k, N_k) {
  (N_k + alpha[k])/(length(z_n.minus.i) + 1 + sum(alpha)-1)
}

precalculate.S.for.all.clusters <- function(x, z_n) {
  num.all.points <- length(z_n)
  num.clusters <- max(z_n)
  res <- list()  
  #for (i in 1:num.clusters) {
  for (i in as.integer(names(table(z_n)))) {
    S<- matrix(c(0,0,0,0),nrow=2)
    for (j in 1:num.all.points) {
      S <- S + as.numeric(z_n==i)[j]*t(matrix(x[j,], nrow=1))%*%(matrix(x[j,], nrow=1))
    }
    res[[i]] <- S
  }
  res
}

remove.cluster.point.from.S <- function(precalculated.S, pnt, point.position, z_n) {
  cluster.num <- z_n[point.position]
  sum.to.update <- precalculated.S[[cluster.num]]
  sum.to.update <- sum.to.update - t(matrix(pnt, nrow=1))%*%(matrix(pnt, nrow=1))
  precalculated.S[[cluster.num]] <- sum.to.update
  precalculated.S
} 

add.cluster.point.to.S <- function(precalculated.S, pnt, cluster.num) {
  sum.to.update <- c()

  if (length(precalculated.S) < cluster.num || is.null(precalculated.S[[cluster.num]])) {
    sum.to.update <- matrix(c(0,0,0,0),nrow=2)
  } else {
    sum.to.update <- precalculated.S[[cluster.num]]
  }
  sum.to.update <- sum.to.update + t(matrix(pnt, nrow=1))%*%(matrix(pnt, nrow=1))
  precalculated.S[[cluster.num]] <- sum.to.update
  precalculated.S
}


precalculate.point.sum.for.all.clusters <- function(x, z_n) {
  num.all.points <- length(z_n)
  num.clusters <- max(z_n)
  res <- list()  
  #for (i in 1:num.clusters) {
  for (i in as.integer(names(table(z_n)))) {
    cur.sum <- c(0,0)
    for (j in 1:num.all.points) {
      cur.sum <- cur.sum + as.numeric(z_n==i)[j]*x[j,]
    }
    res[[i]] <- cur.sum
  }
  res
}

remove.cluster.point.from.point.sum <- function(precalculated.point.sum, pnt, point.position, z_n) {
  cluster.num <- z_n[point.position]
  sum.to.update <- precalculated.point.sum[[cluster.num]]
  sum.to.update <- sum.to.update - pnt
  precalculated.point.sum[[cluster.num]] <- sum.to.update
  precalculated.point.sum
} 


add.cluster.point.to.point.sum <- function(precalculated.point.sum, pnt, cluster.num) {
  sum.to.update <- c()
  if (length(precalculated.point.sum) < cluster.num || is.null(precalculated.point.sum[[cluster.num]])) {
    sum.to.update <- c(0,0)
  } else {
    sum.to.update <- precalculated.point.sum[[cluster.num]]
  }
  sum.to.update <- sum.to.update + pnt
  precalculated.point.sum[[cluster.num]] <- sum.to.update
  precalculated.point.sum
}

#http://stats.stackexchange.com/questions/68476/drawing-from-the-multivariate-students-t-distribution

run.collapsed.gibbs <- function(res, num.clusters, num.iters, progress.fun = NULL) {
  
  # Init data
  x <- res$generated.data$data.points
  #number of points   
  n <- res$input.params$num.data.points
  
  m_0 <- res$input.params$m_0
  nu_0 <- res$input.params$nu_0
  V_0 <-  res$input.params$V_0
  S_0 <- solve(res$input.params$S_0.inv)
  alpha_k <- res$input.params$alpha
  
  # For experiments
  alpha_k <- rep(1, num.clusters)
  
  # for plotting
  plot.list <-list()  
  x.df <- as.data.frame(x)
  colnames(x.df) <- c("x.col","y.col")
  
  # choose an initial z (allocation points to clusters)
  p.out <- rep (1.0/num.clusters, num.clusters)
  z_n <- sample(1:num.clusters, n, prob = p.out, replace = TRUE) 
  
  # the actual allocation that we will try to reproduce
  z_n.data <- res$generated.data$points.to.clusters
  actual.cluster.centers.df <- as.data.frame(res$generated.data$cluster.centers)
  colnames(actual.cluster.centers.df) <- c("x.col","y.col")
  
  for (iter in 1:num.iters) {
    
    cat("######### iter ", iter, "\n")
    precalculated.S <- precalculate.S.for.all.clusters(x, z_n)
    precalculated.point.sum <- precalculate.point.sum.for.all.clusters(x, z_n)
    cur.cluster.point.nums <- table(z_n)
    # Plot for the current iteration
    current.cluster.centers.df <- as.data.frame(t(sweep(t(do.call(rbind, precalculated.point.sum)), 2, cur.cluster.point.nums, "/")))
    colnames(current.cluster.centers.df) <- c("x.col","y.col")
    cur.plot <- plot.cur.results(x.df, z_n, z_n.data, current.cluster.centers.df, actual.cluster.centers.df, "Collapsed Gibbs", iter)
    plot.list[[length(plot.list)+1]] <- cur.plot
    if (iter %% 5 ==0 || iter == 1) {
      grid.arrange(cur.plot, ncol=1)
    }
  
    # n is the number of points
    #permute the points for every iteration
    points.permuted <- sample.int(n,n)
    for (pnt in 1:n) {
      if (pnt %% 100 == 0) {
        cat("point ", pnt, "")
      }
      i <- points.permuted[pnt]  
      
      
      z_n.minus.i <- z_n[-i]
      cluster.probs <- c()
      x_i <- x[i,]
      x.minus.i <- x[-i,,drop=FALSE]
      precalculated.S <- remove.cluster.point.from.S(precalculated.S, x_i, i, z_n)
      precalculated.point.sum <- remove.cluster.point.from.point.sum(precalculated.point.sum, x_i, i, z_n)
      N_k.table <- table(z_n.minus.i)
      for (k in 1:num.clusters) {
        # a cluster may actually disappear completely so we need to be careful
        N_k <- N_k.table[as.character(k)]
        if (!is.na(N_k)) { 
          P_z_k <-get.P_z_k(alpha_k, z_n.minus.i, k, N_k)
          
          #FIX ME: - what is the right param?
          kappa_0 <- 0.2
          #FIX ME: This is the number of dimensions
          D <-2
          
          m_k.S_k_i <- get.m_k.S_k_i(S_0, m_0, kappa_0, nu_0, z_n.minus.i, x.minus.i, k, precalculated.point.sum, precalculated.S, N_k)
          
          kappa_N <- kappa_0 + N_k
          nu_N <- nu_0 + N_k
          # equation (15)
          cluster.probs <- c(cluster.probs, P_z_k*dmt(x_i, mean=m_k.S_k_i[[1]], kappa_N/(kappa_N*(nu_N-D+1)) * m_k.S_k_i[[2]], nu_N-D+1))
        } else {
          cluster.probs <- c(cluster.probs, 0)
        }
      }
      
      sml <- rmultinom(1, size = 1, prob = cluster.probs)
      sml <- sml[,1]
      assigned.cluster <- which.max(sml)
      z_n[i] <- assigned.cluster 
      precalculated.S <- add.cluster.point.to.S(precalculated.S, x_i, assigned.cluster)
      precalculated.point.sum <- add.cluster.point.to.point.sum(precalculated.point.sum, x_i, assigned.cluster)
    }
    if (!is.null(progress.fun)) {
      progress.fun(1/num.iters, detail = paste("Doing iter", iter))
    }
  }
  
  # Plot the final result
  precalculated.point.sum <- precalculate.point.sum.for.all.clusters(x, z_n)
  cur.cluster.point.nums <- table(z_n)
  current.cluster.centers.df <- as.data.frame(t(sweep(t(do.call(rbind, precalculated.point.sum)), 2, cur.cluster.point.nums, "/")))
  colnames(current.cluster.centers.df) <- c("x.col","y.col")
  cur.plot <- plot.cur.results(x.df, z_n, z_n.data, current.cluster.centers.df, actual.cluster.centers.df, "Collapsed Gibbs - final", NULL)
  plot.list[[length(plot.list)+1]] <- cur.plot
  
  list(z_n = z_n, plot.list=plot.list)
}

#z_n.collapsed.gibbs <- run.collapsed.gibbs(res, res$input.params$num.clusters, 30)

#### Infinite Mixture ####

get.P_z_k.infinite.mixture <- function(alpha, z_n.minus.i, N_k) {
  (N_k + alpha)/(length(z_n.minus.i) + 1 + alpha-1)
}

#mew cluster 
get.m.S.new.cluster <- function(S_0, m_0, kappa_0, nu_0, x) {
  
  x_bar_i <- x
  N_i <- 1
  # equation (7)
  m_k_i <- (kappa_0*m_0+N_i*x_bar_i)/(kappa_0+N_i) 
  
  # equation (8)
  # This is 2 by 2 matrix
  S<- matrix(c(0,0,0,0),nrow=2)
  S <- S + t(matrix(x, nrow=1))%*%(matrix(x, nrow=1))
  S_k_i <- S_0 + S 
  S_k_i <- S_k_i + kappa_0*t(matrix(m_0, nrow=1))%*%(matrix(m_0, nrow=1))
  S_k_i <- S_k_i - (kappa_0 + N_i)*t(matrix(m_k_i, nrow=1))%*%(matrix(m_k_i, nrow=1))
  res <- list()
  res[[1]] <- m_k_i
  res[[2]] <- S_k_i
  res
}

run.infinite.mixture.collapsed.gibbs <- function(res, num.iters, progress.fun = NULL) {
  
  # Init data
  #number of points   
  n <- res$input.params$num.data.points
  
  x <- res$generated.data$data.points
  m_0 <- res$input.params$m_0
  nu_0 <- res$input.params$nu_0
  V_0 <-  res$input.params$V_0
  S_0 <- solve(res$input.params$S_0)
  
  #alpha_k <- res$alpha
  alpha <- 0.03  
  
  # For plotting
  plot.list <-list()
  x.df <- as.data.frame(x)
  colnames(x.df) <- c("x.col","y.col")
  
  # the actual allocation that we will try to reproduce 
  z_n.data <- res$generated.data$points.to.clusters
  actual.cluster.centers.df <- as.data.frame(res$generated.data$cluster.centers)
  colnames(actual.cluster.centers.df) <- c("x.col","y.col")
  
  
  init.number.of.clusters <- 1
  z_n <- sample(init.number.of.clusters, n, replace=TRUE)
  
  num.clusters <- init.number.of.clusters
  points.in.clusters <- table(z_n)  
  #since clusters can disappear we will keep the ids of existing clusters
  cluster.ids <- seq(1,num.clusters)
  
  #Set up colors
  #http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
  #http://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
  colorCount <- 22
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  #avail.colors <- sample(rainbow(colorCount))
  avail.colors <- getPalette(colorCount)
  cluster.colors <- avail.colors[1:num.clusters]
  avail.colors <- avail.colors[(num.clusters+1):length(avail.colors)]
  
  
  for (iter in 1:num.iters) {
    
    cat("######### iter ", iter, "\n")
    
    precalculated.S <- precalculate.S.for.all.clusters(x, z_n)
    precalculated.point.sum <- precalculate.point.sum.for.all.clusters(x, z_n)
    
    cur.cluster.point.nums <- table(z_n)
    current.cluster.centers.df <- as.data.frame(t(sweep(t(do.call(rbind, precalculated.point.sum)), 2, cur.cluster.point.nums, "/")))
    colnames(current.cluster.centers.df) <- c("x.col","y.col")  
    # Plot for the current iteration
    cur.plot <- plot.cur.results(x.df, z_n, z_n.data, current.cluster.centers.df, actual.cluster.centers.df, "Infinite Mixture", iter, cluster.colors)
    plot.list[[length(plot.list)+1]] <- cur.plot
    if (iter %% 5 ==0 || iter == 1) {
      grid.arrange(cur.plot, ncol=1)
    }
    print(table(z_n))
    
    
    #permute the points for every iteration
    points.permuted <- sample.int(n,n)  
    for (pnt in 1:n) {
      if (pnt %% 100 == 0) {
        cat("point ", pnt, "")
      }
      
      
      i <- points.permuted[pnt]
      
      # remove point from the cluster allocation list
      cluster.of.point <- z_n[i]
      # adjust the number of points in this cluster
      pos <- which(cluster.ids == cluster.of.point)
      # should be exactly one
      stopifnot(length(pos) == 1)
      if (points.in.clusters[pos] > 1) {
        points.in.clusters[pos] <- points.in.clusters[pos] - 1
      } else {
        # we need to remove this cluster
        cat(" need to remove a cluster ")
        num.clusters <- num.clusters-1
        points.in.clusters <- points.in.clusters[-pos]
        cluster.ids <- cluster.ids[-pos]
        cluster.color <- cluster.colors[pos]
        cluster.colors <- cluster.colors[-pos]
        avail.colors <- c(avail.colors, cluster.color)
      }
      
      x_i <- x[i,]
      precalculated.S <- remove.cluster.point.from.S(precalculated.S, x_i, i, z_n)
      precalculated.point.sum <- remove.cluster.point.from.point.sum(precalculated.point.sum, x_i, i, z_n)
      
      z_n.minus.i <- z_n[-i]
      cluster.probs <-c()
      N_k.table <- table(z_n.minus.i) 
      for (cluster.pos in 1:length(cluster.ids)) {
        k <- cluster.ids[cluster.pos]
        N_k <- N_k.table[cluster.pos]
        P_z_k <-get.P_z_k.infinite.mixture(alpha, z_n.minus.i, N_k)
        x_i <- x[i,]
        x.minus.i <- x[-i,,drop=FALSE]
        
        # Better way to set kappa_0?
        kappa_0 <- 0.2
        # dimensions
        D <-2
        
        m_k.S_k_i <- get.m_k.S_k_i(S_0, m_0, kappa_0, nu_0, z_n.minus.i, x.minus.i, k, precalculated.point.sum, precalculated.S, N_k)
        
        kappa_N <- kappa_0 + N_k
        nu_N <- nu_0 + N_k
        # equation (15) + line 9 of Algorithm 1 in Kamper
        cluster.probs <- c(cluster.probs, P_z_k*dmt(x_i, mean=m_k.S_k_i[[1]], kappa_N/(kappa_N*(nu_N-D+1)) * m_k.S_k_i[[2]], nu_N-D+1))  
      }
      # add probs for a new cluster
      P_z_k.new.cluster <- alpha/(length(z_n.minus.i) + 1 + alpha-1)
      m_k.S_k.new.cluster <- get.m.S.new.cluster(S_0, m_0, kappa_0, nu_0, x_i)
      kappa_N.new.cluster <- kappa_0 + 1
      nu_N.new.cluster <- nu_0 + 1
      cluster.probs <- c(cluster.probs, P_z_k.new.cluster*dmt(x_i, mean=m_k.S_k.new.cluster[[1]], kappa_N.new.cluster/(kappa_N.new.cluster*(nu_N.new.cluster-D+1)) * m_k.S_k.new.cluster[[2]], nu_N.new.cluster-D+1))  
      sml <- rmultinom(1, size = 1, prob = cluster.probs)
      sml <- sml[,1]
      chosen.cluster.pos <- which.max(sml)
      if (chosen.cluster.pos == length(cluster.probs)) {
        #new cluster was selected
        cat(" new cluster created ")
        num.clusters <- num.clusters+1
        cluster.ids <- c(cluster.ids, max(cluster.ids)+1)
        points.in.clusters<-c(points.in.clusters,1)
        cluster.colors <- c(cluster.colors, avail.colors[1])
        avail.colors <- avail.colors[-1]
      } 
      assigned.cluster <- cluster.ids[chosen.cluster.pos] 
      z_n[i] <- assigned.cluster 
      precalculated.S <- add.cluster.point.to.S(precalculated.S, x_i, assigned.cluster)
      precalculated.point.sum <- add.cluster.point.to.point.sum(precalculated.point.sum, x_i, assigned.cluster)
      points.in.clusters <- table(z_n)
    }
    if (!is.null(progress.fun)) {
      progress.fun(1/num.iters, detail = paste("Doing iter", iter))
    }
  }
  
  # Plot the final result
  precalculated.point.sum <- precalculate.point.sum.for.all.clusters(x, z_n)
  cur.cluster.point.nums <- table(z_n)
  current.cluster.centers.df <- as.data.frame(t(sweep(t(do.call(rbind, precalculated.point.sum)), 2, cur.cluster.point.nums, "/")))
  colnames(current.cluster.centers.df) <- c("x.col","y.col")  
  # Plot for the current iteration
  cur.plot <- plot.cur.results(x.df, z_n, z_n.data, current.cluster.centers.df, actual.cluster.centers.df, "Infinite Mixture - final", NULL, cluster.colors)
  plot.list[[length(plot.list)+1]] <- cur.plot
  
  
  list(z_n = z_n, plot.list=plot.list)
}

### Profiling
#rprof.filename <-"infinite.mixture.rprof"
#Rprof(rprof.filename)
#z_n.infinite.mixture <- run.infinite.mixture.collapsed.gibbs(res, 60)
#Rprof(NULL)

#rprof.summary <- summaryRprof(rprof.filename)
#rprof.summary$by.self
#rprof.summary$by.total


#### K-Means ####

#http://r.789695.n4.nabble.com/k-means-td4274383.html
require(pracma)

kmpp <- function(X, k) {
  n <- nrow(X)
  centroids <- numeric(k)
  centroids[1] <- sample(1:n, 1)
  
  if (k > 1) {
    for (i in 2:k) {
      dm <- distmat(X, X[centroids, ])
      pr <- apply(dm, 1, min); pr[centroids] <- 0
      centroids[i] <- sample(1:n, 1, prob = pr)
    }
  }
  X[centroids, ,drop=FALSE]
} 

kmeans.cluster.points <- function(data.points, centroids) {
  dm <- distmat(data.points, centroids)
  cluster.assignment <- apply(dm, 1, which.min)
  cluster.assignment
}

kmeans.recompute.centroids <- function(data.points, cluster.assignments, num.clusters) {
  new.centroids <- c()
  for(i in 1:num.clusters) {
    new.centroid <- colMeans(data.points[cluster.assignments == i,,drop=FALSE])
    new.centroids <- rbind(new.centroids, new.centroid)
  } 
  rownames(new.centroids) <- NULL
  new.centroids
}

run.kmeans <- function(res, num.clusters, num.iterations, progress.fun = NULL) {
  print("Running k-means")
  x <- res$generated.data$data.points
  
  # For plotting
  plot.list <-list()
  x.df <- as.data.frame(x)
  colnames(x.df) <- c("x.col","y.col")
  actual.cluster.centers.df <- as.data.frame(res$generated.data$cluster.centers)
  colnames(actual.cluster.centers.df) <- c("x.col","y.col")
  z_n.data <- res$generated.data$points.to.clusters  
  
  z_n <- numeric(size(x,1))
  # find the initial centroids
  init.centroids <- kmpp(x, num.clusters)
  
  centroids <- init.centroids
  
  for (iter in 1:num.iterations) { 
    
    # Generate iteration plot
    centroids.df <- as.data.frame(centroids)
    colnames(centroids.df) <- c("x.col","y.col")                              
    cur.plot <- plot.cur.results(x.df, z_n, z_n.data, centroids.df, actual.cluster.centers.df, "K-Means", iter)
    plot.list[[length(plot.list)+1]] <- cur.plot 
    if (iter %% 2 == 0 || iter == 1) {
      grid.arrange(cur.plot, ncol=1)
    }
    
    z_n.prev <- z_n
    z_n <- kmeans.cluster.points(x, centroids) 
    if (sum(abs(z_n-z_n.prev)) == 0) {
      cat("No change - aborting iterations")
      break
    }
    centroids <- kmeans.recompute.centroids(x, z_n, num.clusters)
    
    if (!is.null(progress.fun)) {
      progress.fun(1/num.iterations, detail = paste("Doing iter", iter))
    }
    
  }
  # Generate the final result plot
  centroids.df <- as.data.frame(centroids)
  colnames(centroids.df) <- c("x.col","y.col")                              
  cur.plot <- plot.cur.results(x.df, z_n, z_n.data, centroids.df, actual.cluster.centers.df, "K-Means - final", NULL)
  plot.list[[length(plot.list)+1]] <- cur.plot 
  
  list(z_n = z_n, plot.list=plot.list)
}

# test #

#centroids <- kmpp(res$generated.data$data.points, res$input.params$num.clusters)
#as.data.frame(kmpp(res$generated.data$data.points, 1))
#dm <- distmat(res$generated.data$data.points, centroids)
#cluster.assignment <- apply(dm, 1, which.min)

# end test #

#run.kmeans(res, res$input.params$num.clusters, 10)
#res.kmeans <- run.kmeans(res, 5, 10)
#generate.animation(res.kmeans)


### Subspace clustering ####
# LAC (Locally Adaptive Clustering) 
# See "Subspace Clustering of High Dimensional Data" by Carlotta Domeniconi, Dimitris Papadopoulos, Dimitrios Gunopulos, Sheng Ma. 

lac.cluster.points <- function(data.points, centroids, weights) {
  num.points <- size(data.points,1)
  cluster.assignment <- numeric(num.points)
  num.clusters <- size(centroids, 1)
  for (i in 1:num.points) {
    rep.pnt <- matrix(rep(data.points[i,], num.clusters), nrow=num.clusters, byrow=TRUE)
    tmp <- (rep.pnt-centroids)^2*weights
    dst.from.centroids <- rowSums(tmp)
    cluster.assignment[i] <- which.min(dst.from.centroids)
  }
  cluster.assignment
}

lac.recompute.weights <- function(data.points, centroids, z_n, h, N) {
  num.points.in.clusters <- table(z_n)
  num.clusters <- size(centroids, 1)
  num.points <- size(data.points, 1)
  X_ji <-matrix(rep(0, num.clusters*N), nrow=num.clusters)
  for (j in 1:num.clusters) {
    # cluster my disappear
    if (num.points.in.clusters[j] == 0) {
      # this will make the weights huge
      X_ji[j,] = rep(-Inf, N)
      next
    }
    for (i in 1:N) {
      X_ji[j,i] <-0
      for (k in 1:num.points) {
        X_ji[j,i] = X_ji[j,i]  + as.numeric(z_n==j)[k]*(data.points[k,i]-centroids[j,i])^2
      }
      X_ji[j,i] <- X_ji[j,i]/num.points.in.clusters[j]
    }
  }
  w_ji <-matrix(rep(0, num.clusters*N), nrow=num.clusters)
  for (j in 1:num.clusters) {
    for (i in 1:N) {
      w_ji[j,i] <- exp(-h*X_ji[j,i]) 
    }
  }
  
  for (j in 1:num.clusters) {
    X_jl <- 0
    for (l in 1:N) {
      X_jl <- X_jl + exp(-h*2*X_ji[j,l])
    }
    w_ji[j,] <- w_ji[j,]/sqrt(X_jl)
  }
  
  w_ji
}

lac.recompute.centroids <- function(data.points, cluster.assignments, num.clusters) {
  new.centroids <- c()
  for(i in 1:num.clusters) {
    new.centroid <- colMeans(data.points[cluster.assignments == i,,drop=FALSE])
    new.centroids <- rbind(new.centroids, new.centroid)
  } 
  rownames(new.centroids) <- NULL
  new.centroids
}


run.lac <- function(res, num.clusters, num.iterations, progress.fun = NULL) {
  print("Running lac")
  x <- res$generated.data$data.points
  # for plotting
  x.df <- as.data.frame(x)
  colnames(x.df) <- c("x.col","y.col")
  
  # For plotting
  plot.list <-list()
  actual.cluster.centers.df <- as.data.frame(res$generated.data$cluster.centers)
  colnames(actual.cluster.centers.df) <- c("x.col","y.col")
  z_n.data <- res$generated.data$points.to.clusters  
  
  z_n <- numeric(size(x,1))
  # find the initial centroids
  init.centroids <- kmpp(x, num.clusters)
  
  # initialize weights
  # number of dimensions
  N <- 2
  h <- 2
  w_ji <- matrix(rep(sqrt(N), num.clusters*N), nrow=num.clusters)
  
  centroids <- init.centroids
  
  for (iter in 1:num.iterations) { 
    
    # Plot for the current iteration
    centroids.df <- as.data.frame(centroids)
    colnames(centroids.df) <- c("x.col","y.col")                              
    cur.plot <- plot.cur.results(x.df, z_n, z_n.data, centroids.df, actual.cluster.centers.df, "LAC", iter)
    plot.list[[length(plot.list)+1]] <- cur.plot
    
    if (iter %% 2 == 0 || iter == 1) {
      grid.arrange(cur.plot, ncol=1)
    }
    
    z_n.prev <- z_n
    z_n <- lac.cluster.points(x, centroids, w_ji) 
    if (sum(abs(z_n-z_n.prev)) == 0) {
      cat("No change - aborting iterations")
      break
    }
    w_ji <- lac.recompute.weights(x, centroids, z_n, h, N)
    z_n <- lac.cluster.points(x, centroids, w_ji)
    centroids <- lac.recompute.centroids(x, z_n, num.clusters)
    
    if (!is.null(progress.fun)) {
      progress.fun(1/num.iterations, detail = paste("Doing iter", iter))
    }
  }
  
  # Plot the final result
  centroids.df <- as.data.frame(centroids)
  colnames(centroids.df) <- c("x.col","y.col")                              
  cur.plot <- plot.cur.results(x.df, z_n, z_n.data, centroids.df, actual.cluster.centers.df, "LAC - final", NULL)
  plot.list[[length(plot.list)+1]] <- cur.plot
  
  
  list(z_n = z_n, plot.list=plot.list)
}


#num.clusters <- 9
#num.data.points <- 80
#boundary <-20
#res <- generate.data(num.clusters, num.data.points, boundary, TRUE)
#plot.generated.data(res)
#run.lac(res, res$input.params$num.clusters, 20)

# setwd("<file-loc-here>")
# runApp(host = host.name, port=5533)

