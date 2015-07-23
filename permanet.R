# create some data

set.seed(10)
n_nodes <- 100
nodes <- data.frame(id = 1:n_nodes,
                    sex = sample(c("F", "M"), n_nodes, replace=T),
                    grid = sample(LETTERS[1:2], n_nodes, replace=T))

# model for network
p_sex <- matrix(c(0.5,1,1,0.5), 2, 2) # twice as likely contact for dyad
p_grid <- matrix(c(1,0,0,1), 2, 2) # no cross-grid
p_contact <- 0.3

# model for contacts on network
c_sex <- matrix(c(1,2,2,1), 2, 2) # twice as likely number of contacts for dyad
c_grid <- matrix(c(2,0,0,3), 2, 2) # grid B more likely contact than grid A

# simulate network
network <- matrix(NA, nrow(nodes), nrow(nodes))
rownames(network) <- nodes$id
colnames(network) <- nodes$id
for (i in nodes$id) {
  for (j in nodes$id) {
    p <- p_sex[nodes$sex[i], nodes$sex[j]] * p_grid[nodes$grid[i], nodes$grid[j]] * p_contact
    log_c <- c_sex[nodes$sex[i], nodes$sex[j]] + c_grid[nodes$grid[i], nodes$grid[j]]
    network[i,j] <- rpois(1, p*exp(log_c))
  }
}
# make symmetric
for (i in seq_len(nrow(network)-1)) {
  for (j in i + seq_len(nrow(network)-i)) {
    network[i,j] <- network[j,i] <- network[i,j] + network[j,i]
  }
  network[i,i] <- 0
}

# righto. We want to see if the contacts are dependent on sex, in particular if dyad
# is different to same-sex, accounting for grid.

# thus, we work within-grid to work out the average contact of same-sex vs different-sex
test_stat <- function(network, nodes) {

  # Note: this assumes that net and nodes are indexed similarly
  test_stat_inner <- function(net, nodes) {
    diff <- 0
    for (i in 1:nrow(net)) {
      for (j in 1:nrow(net)) {
        contact <- net[i,j]
        diff <- diff + contact * ((nodes$sex[i] != nodes$sex[j]) - (nodes$sex[i] == nodes$sex[j]))
      }
    }
    diff
  }

  outer <- tapply(1:nrow(nodes), nodes$grid, c)
  stat <- 0
  for (i in 1:length(outer)) {
    stat <- stat + test_stat_inner(network[outer[[i]], outer[[i]]], nodes[outer[[i]],])
  }
  stat
}

# thus, we work within-grid to work out the average contact of same-sex vs different-sex
test_stat2 <- function(network, group, cov_levels, cov_map) {
  
  # Note: this assumes that net and nodes are indexed similarly
  test_stat_inner <- function(net, cov_levels, cov_map) {
    sums <- numeric(length(cov_levels))
    for (i in 1:length(cov_levels)) {
      sums[i] <- sum(net[cov_map == cov_levels[i]])
    }
    sums[2] - sums[1]
    # f-statistic will be different, we'll want:
    # 1. variance within each group
    # 2. mean within each group
    # 3. 
  }

  # Note: this assumes that net and nodes are indexed similarly
  test_stat_inner <- function(net, cov_levels, cov_map) {
    sums <- numeric(length(cov_levels))
    for (i in 1:length(cov_levels)) {
      sums[i] <- sum(net[cov_map == cov_levels[i]])
    }
    sums[2] - sums[1]
    # f-statistic will be different, we'll want:
    # 1. variance within each group
    # 2. mean within each group
    # 3. 
  }
  
  outer <- tapply(1:nrow(nodes), group, c)
  stat <- 0
  for (i in 1:length(outer)) {
    stat <- stat + test_stat_inner(network[ outer[[i]], outer[[i]] ],
                                   cov_levels,
                                   cov_map[ outer[[i]], outer[[i]] ])
  }
  stat
}

# thus, we work within-grid to work out the average contact of same-sex vs different-sex
test_stat3 <- function(network, group, cov_levels, cov_map) {

  test_stat_inner <- function(net, cov_levels, cov_map) {
    #    sums <- numeric(length(cov_levels))
    g <- length(cov_levels)
    means <- numeric(g)
    vars  <- numeric(g)
    n     <- numeric(g)
    for (i in 1:g) {
      means[i] <- mean(net[cov_map == cov_levels[i]])
      vars[i] <- var(as.numeric(net[cov_map == cov_levels[i]]))
      n[i] <- sum(cov_map == cov_levels[i]) # NOTE: This will be constant
    }
    var_g <- var(means)
    var_r <- sum((n-1)*vars)/sum((n-1)) # weighted avg of residual variance
    (var_g / (g-1)) / (var_r / (sum(n)-g))
  }
  
  outer <- tapply(1:nrow(nodes), group, c)
  stat <- 0
  for (i in 1:length(outer)) {
    stat <- stat + test_stat_inner(network[ outer[[i]], outer[[i]] ],
                                   cov_levels,
                                   cov_map[ outer[[i]], outer[[i]] ])
  }
  stat
}

e <- expand.grid(nodes$sex, nodes$sex)
sex_net <- matrix(ifelse(e[,1] == e[,2],1,2), nrow(nodes), nrow(nodes))
sex_lev <- 1:2

ts_data <- test_stat2(network, nodes$grid, sex_lev, sex_net)

# now we randomly permute the possums within the grid
num_iter <-1000
outer <- tapply(1:nrow(nodes), nodes$grid, c)
ord <- matrix(NA, num_iter, nrow(nodes))
for (it in seq_len(num_iter)) {
  for (i in 1:length(outer)) {
    ord[it,outer[[i]]] <- (1:nrow(nodes))[sample(outer[[i]])]
  }
}

ts <- numeric(num_iter)
for (it in seq_len(num_iter)) {
  o <- ord[it,]
  ts[it] <- test_stat(network, nodes[o,])
}

ts2 <- numeric(num_iter)
for (it in seq_len(num_iter)) {
  o <- ord[it,]
  ts2[it] <- test_stat2(network, nodes$grid[o], sex_lev, sex_net[o,o])
}

any(ts != ts2)

ts3_data <- test_stat3(network, nodes$grid, sex_lev, sex_net)
ts3 <- numeric(num_iter)
for (it in seq_len(num_iter)) {
  o <- ord[it,]
  ts3[it] <- test_stat3(network, nodes$grid[o], sex_lev, sex_net[o,o])
}

# do a histogram
hist(ts3)
abline(v=ts_data, col="red")
