# Calculate power analyses for common language effect size
setwd("~/Documents/Projects/Ammonites/code/")


library(auRoc)
library(ggplot2)
library(reshape2)
library(dplyr)
library(forcats)

# Specify plotting aesthetics
darkblue  <- '#12255A'
orange    <- '#E38C59'
errbarcol <- '#CA282C'
font <- 'Palatino'
fsize <- 9


# This simulation assumes data is sampled from two independent normal 
# distributions. The number of samples and distribution variances are taken from 
# the dataset for which the power analysis is being conducted. 

# Effect sizes, expressed in terms of common language effect sizes, also known 
# as the probability of superiority, are systematically explored from 0.5+ to 1.
# This is done by shifting one of the normal distributions along the horizontal
# axis, changing the fraction that has larger values. 

###################  Functions ################### 

# Take an effect size and calculate what the mean of the shifting distribution
# needs to be

find_mean_shift <- function(theta, x.mean, x.var, y.var) {
  # Find how by how much the mean of the normal distribution Y has to be 
  # shifted so that normal distribution Y has a probability of superiority of 
  # theta relative to X. 
  # Variables:
  #   theta   = Common language effect size / probability of superiority
  #   x.mean  = Mean of the non-shifting distribution X
  #   x.var   = Variation of non-shifting distribution X
  #   y.var   = Variation of shifting distribution Y
  # Value:
  #   y.mean  = The mean of shifting distribution Y 
  
  # Calculate the difference between the two normal distributions.
  # This difference is itself normally distributed. Theta specifies the probability
  # that the difference variable exceeds zero (i.e. the values in normal
  # distribution X are larger than values in normal distribution Y).
  
  # Assuming that the variances are uncorrelated, they can simply be added
  diff.var = x.var + y.var
  
  # Use the quantile function of a standard normal distribution, convert by
  # multiplying quantile function by standard deviation of difference distribution
  # I am looking for _larger than_, so lower.tail = FALSE
  # It is already known that x, the value to be exceeded, is zero
  # x = mu + sigma * qnorm( theta )
  diff.mean <- sqrt(diff.var) * qnorm(theta, mean = 0, sd = 1, lower.tail = FALSE)
  y.mean <- x.mean - diff.mean
  return(y.mean)
}


determine_power <- function(theta, m, n, x.var, y.var, n.samples = 10000,
  ci.method = 'DL.corr', alpha.level = 0.05, nboot = 1000) {
  # Simulate two distributions that are shifted relatively to one another, so 
  # that distribution Y has a probability of superiority of theta over
  # distribtion X. For a total of n.samples, take m and n random samples from
  # both distributions. The distributions have variances x.var and y.var 
  # respectively.
  # Variables:
  #   theta     = Common language effect size / probability of superiority
  #   m         = Number of samples taken from distribution X
  #   n         = Number of samples taken from distribution Y
  #   x.var     = Variation of non-shifting distribution X
  #   y.var     = Variation of shifting distribution Y
  #   n.samples = Number of sampling replicates
  #   ci.method = Confidence interval method from auRoc. Must be one of 
  #               "newcombe", "pepe", "delong", "jackknife", "bootstrapP",
  #               "bootstrapBCa"
  #   alpha     = Alpha level. Default 0.05
  #   nboot     = Number of bootstrap repetitions for auRoc function
  # Value:
  #   power     = Fraction of true positives, or 1 - rate of false negatives
  
  x.mean <- 0 # This distribution does not shift. Exact value of mean irrelevant

  # Find how much the mean of distribution y has to shift to get a probability
  # of superiority equivalent to theta
  y.mean <- find_mean_shift(theta, x.mean, x.var, y.var)

  
  # Sample 10,000 trials with sizes m & n from the distributions
  x.samples <- rnorm(m*n.samples, x.mean, sqrt(x.var))
  y.samples <- rnorm(n*n.samples, y.mean, sqrt(y.var))
  # Reshape
  x.trials <- matrix(x.samples, nrow = n.samples, byrow = TRUE)
  y.trials <- matrix(y.samples, nrow = n.samples, byrow = TRUE)
  
  # Is a significant difference detected or not?
  detected <- c()
  effsize <- c()
  
  for (i in 1:n.samples) {
    # Calculate effect size and confidence intervals
    e.size <- auc.nonpara.mw(y.trials[i,], x.trials[i,], 
      method = ci.method, conf.level = 1 - alpha.level, nboot = nboot)
    
    # Determine if the confidence interval crosses the 0.5 mark. If it does, no 
    # significant effect detected. If not, mark as significant. 
    # Significant: Either both are smaller, or both are larger than 0.5
    if ((e.size[2] < 0.5 & e.size[3] < 0.5) | (e.size[2] > 0.5 & e.size[3] > 0.5)) 
    {
      detected <- append(detected, TRUE)
    } else {
      detected <- append(detected, FALSE)
    }
    effsize <- append(effsize, e.size[1]) # Sanity check! Effect sizes == theta
  }
  # Power is 1 - rate of false negatives
  power <- sum(detected) / n.samples
  return(power)
}

###################  Run analysis ################### 

# Read in empirical data
effect.sizes.ammon <- 
  read.csv('../results/comparing_hypotheses/ammonoids_effect_sizes.csv')
effect.sizes.nauti <-
  read.csv('../results/comparing_hypotheses/nautilids_effect_sizes.csv')

# Standardize variances to be relative (so that means become irrelevant and 
# the non-shifting normal distribution X can have a mean of zero)
effect.sizes.ammon$std.var.ext <- 
  effect.sizes.ammon$var.ext / effect.sizes.ammon$var.ext
effect.sizes.ammon$std.var.surv <-
  effect.sizes.ammon$var.surv / effect.sizes.ammon$var.ext

effect.sizes.nauti$std.var.ext <- 
  effect.sizes.nauti$var.ext / effect.sizes.nauti$var.ext
effect.sizes.nauti$std.var.surv <-
  effect.sizes.nauti$var.surv / effect.sizes.nauti$var.ext

# Calculate power
simef <- seq(0.5, 0.9, by = 0.1) # simulated effect sizes
am.power <- matrix(0, nrow = nrow(effect.sizes.ammon), ncol = length(simef))
na.power <- matrix(0, nrow = nrow(effect.sizes.nauti), ncol = length(simef))
am.idx <- expand.grid(1:nrow(effect.sizes.ammon), 1:length(simef)) # Indicies of matrix
na.idx <- expand.grid(1:nrow(effect.sizes.nauti), 1:length(simef)) # Indicies of matrix

# Ammonoids
# Progress bar... this is a bit slow
pbar <- txtProgressBar(min = 0, max = nrow(am.idx))

for (i in 1:nrow(am.idx)) {
  
  m     <- effect.sizes.ammon$num.ext[am.idx[i,1]]
  n     <- effect.sizes.ammon$num.surv[am.idx[i,1]]
  x.var <- effect.sizes.ammon$std.var.ext[am.idx[i,1]]
  y.var <- effect.sizes.ammon$std.var.surv[am.idx[i,1]]
  
  am.power[am.idx[i,1], am.idx[i,2]] <- determine_power(simef[am.idx[i,2]], m, n, x.var, y.var)
  
  setTxtProgressBar(pbar, i)
}

# Convert to dataframe
am.power.df <- data.frame(am.power)
colnames(am.power.df) <- simef
am.power.df <- cbind(data.frame(variable = effect.sizes.ammon$variable), am.power.df)

# Save power dataframe
write.csv(am.power.df, '../results/comparing_hypotheses/ammonoids_power.csv',
  row.names = FALSE)

# Nautilids
pbar <- txtProgressBar(min = 0, max = nrow(na.idx))

for (i in 1:nrow(na.idx)) {
  
  m     <- effect.sizes.nauti$num.ext[na.idx[i,1]]
  n     <- effect.sizes.nauti$num.surv[na.idx[i,1]]
  x.var <- effect.sizes.nauti$std.var.ext[na.idx[i,1]]
  y.var <- effect.sizes.nauti$std.var.surv[na.idx[i,1]]
  
  na.power[na.idx[i,1], na.idx[i,2]] <- determine_power(simef[na.idx[i,2]], m, n, x.var, y.var)
  
  setTxtProgressBar(pbar, i)
}

# Convert to dataframe
na.power.df <- data.frame(na.power)
colnames(na.power.df) <- simef
na.power.df <- cbind(data.frame(variable = effect.sizes.nauti$variable), na.power.df)

# Save power dataframe
write.csv(na.power.df, '../results/comparing_hypotheses/nautilids_power.csv',
          row.names = FALSE)


###################  Plotting heatmap ################### 

# Plot a heatmap of power vs effect size for all variables of interest

# Ammonoids
am.power.df.plt <- melt(am.power.df)
colnames(am.power.df.plt) <- c('variable', 'eff.size', 'power')

am.power.df.plt$variable <- fct_rev(factor(am.power.df.plt$variable, 
  levels = effect.sizes.ammon$variable))

# Reframe effect sizes so that they are all > 0.5 
plot.ef.ammon <- effect.sizes.ammon %>%
  mutate(rel.eff.size = if_else(eff.size < 0.5, 1 - eff.size, eff.size))
# rounded for plotting
plot.ef.ammon$eff.size <- 
  as.factor(round(plot.ef.ammon$rel.eff.size, digits = 1))

p <- ggplot(am.power.df.plt, aes(x = eff.size, y = variable)) +
  geom_tile(aes(fill = power)) +
  scale_fill_distiller(name = 'Power', palette = 'Reds', direction = 1,
    limits = c(0,1)) +
  geom_point(data = plot.ef.ammon, aes(x = eff.size, y = variable)) +
  labs(x = 'Effect size', y = '') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    legend.title = element_text(size = fsize, colour = 'black', family = font),
    legend.text = element_text(size = fsize, colour = 'black', family = font),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
  )
ggsave(
  '../results/comparing_hypotheses/Power_analysis_ammonoids.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)

# Nautilids
na.power.df.plt <- melt(na.power.df)
colnames(na.power.df.plt) <- c('variable', 'eff.size', 'power')

na.power.df.plt$variable <- fct_rev(factor(na.power.df.plt$variable, 
  levels = effect.sizes.nauti$variable))

# Reframe effect sizes so that they are all > 0.5 
plot.ef.nauti <- effect.sizes.nauti %>%
  mutate(rel.eff.size = if_else(eff.size < 0.5, 1 - eff.size, eff.size))
# rounded for plotting
plot.ef.nauti$eff.size <- 
  as.factor(round(plot.ef.nauti$rel.eff.size, digits = 1))

p <- ggplot(na.power.df.plt, aes(x = eff.size, y = variable)) +
  geom_tile(aes(fill = power)) +
  scale_fill_distiller(name = 'Power', palette = 'Reds', direction = 1, 
    limits = c(0, 1)) +
  geom_point(data = plot.ef.nauti, aes(x = eff.size, y = variable)) +
  labs(x = 'Effect size', y = '') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = fsize, colour = 'black', family = font),
    axis.text.y = element_text(size = fsize, colour = 'black', family = font),
    axis.title.x = element_text(size = fsize, colour = 'black', family = font),
    legend.title = element_text(size = fsize, colour = 'black', family = font),
    legend.text = element_text(size = fsize, colour = 'black', family = font),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
  )
ggsave(
  '../results/comparing_hypotheses/Power_analysis_nautilids.png',
  width = 12, height = 6, units = 'cm', dpi = 600, plot = p)
