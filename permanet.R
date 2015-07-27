# create some simulated data

create_data <- function() {
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
}


####
# read in data for possums
library(dplyr)

node_data <- read.csv("cont_level_for_mod.csv", stringsAsFactors=FALSE)

# FIXME: possum ID 15 is in grid C, so there's some errors in the data
node_data$grid[node_data$pos_e == 15] <- "CD"
node_data$grid[node_data$pos_e == 42 & node_data$id == 34] <- "CD"
node_data$grid[node_data$pos_e == 61 & node_data$id == 113] <- "CD"
# FIXME:

# find the unique possums and their sex and grid
node_test   <- node_data %>% select(id, sex, grid) %>% mutate(grid=ifelse(grid != "CD", grid, NA))
node_second <- node_data %>% select(id=pos_e, sex=sex_e, grid) %>% mutate(grid=ifelse(grid != "CD", grid, NA))

node_test <- rbind(node_test, node_second) %>% unique %>% arrange(id) %>% mutate(sex = ifelse(sex==1, "M", "F"))

node_test <- node_test %>% group_by(id) %>% summarise(sex = ifelse(length(unique(sex)) == 1, sex[1], "NA"), grid = ifelse(length(unique(na.omit(grid))) == 1, unique(na.omit(grid))[1], "NA"))

# finally combine groups C+D as there's not really another suitable way of dealing with them
nodes <- node_test %>% mutate(grid = ifelse(grid == "C" | grid == "D", "CD", grid))

# now the network
network <- matrix(0, nrow(node_test), nrow(node_test))
rownames(network) <- node_test$id
colnames(network) <- node_test$id
for (i in 1:nrow(node_data)) {
  row <- node_data[i,]
  network[as.character(row$id), as.character(row$pos_e)] <- row$contacts
  network[as.character(row$pos_e), as.character(row$id)] <- row$contacts
}
sum(network > 0)
####

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

# network is the network of contacts
# group is the grouping (permutation of observations is done within group)
# effect is the variable under study at the node level (n levels)
# we want to know whether the n^2 - choose(n,2)/2 combinations at the node level differ.
permanet <- function(network, group, effect, iterations = 1000) {

  n_obs <- length(group)

  # create an adjacency network on the effect variable
  # all combinations of effect variable are considered (i.e. 1,2,3 would produce 11,12,13,22,23,33)
  eff_long <- expand.grid(effect, effect)
  eff_long <- factor(apply(eff_long,1,function(x) { paste0(sort(x),collapse="") }))
  eff_net  <- matrix(eff_long, n_obs, n_obs)

  # f-statistic on the data
  f_data <- test_stat3(network, group, levels(eff_long), eff_net)

  # generate some random iterations
  ord <- matrix(NA, iterations, n_obs)
  for (it in seq_len(iterations)) {
  }

  # observations in each group
  group_obs <- tapply(1:n_obs, group, c)
  f_sim  <- numeric(iterations)
  for (it in seq_len(iterations)) {
    # permute the observations within each group
    ord <- numeric(n_obs)
    for (i in 1:length(group_obs)) {
      ord[group_obs[[i]]] <- (1:n_obs)[sample(group_obs[[i]])]
    }
    # and evaluate our f-statistic for the simulation
    f_sim[it] <- test_stat3(network, group[ord], levels(eff_long), eff_net[ord,ord])
    if (it %% 100 == 0) {
      cat("up to iteration", it, "of", iterations, "\n")
    }
  }
  list(p=sum(f_data <= f_sim)/iterations, f_data=f_data, f_sim=f_sim)
}

f_st <- permanet(log(network+1), nodes$grid, nodes$sex)
print(f_st$p)

# do a histogram
hist(f_st$f_sim)
abline(v=f_st$f_data, col="red")

# do a regression on original data
node_data <- node_data %>% mutate(jm_grid = ifelse(grid == "A" | grid == "B", grid, "CD"))

anova(lm(log(contacts+1) ~ jm_grid + as.factor(dyad), data=node_data))
