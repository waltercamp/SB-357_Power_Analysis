# Title: SB-357 Power Analysis
# Author: Walter Campbell
# Initial Date: 07/02/2025
# Last Edited Date: 07/01/2026
# Purpose: To produce an estimate of power for the 
# SB-357g Study, 960000-1701-000-01262.

# Clear workspace
rm(list=ls())

# Install and load packages - More than I ended up using, but can't hurt to have all
packages <- c(
  "arrow", "dataReporter", "data.table", "fs", "Hmisc", "here", "haven",
  "janitor", "lubridate", "magrittr", "purrr", "readr", "skimr", "stringr",
  "tidyr", "tidyselect", "labelled", "dplyr", "DescTools", "broom", 
  "janitor", "gmodels", "crosstable", "flextable", "pwr", "ggplot2", "SimEngine", 
  "tidyverse", "pwrss", "effectsize"
)
install_missing_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) install.packages(new_packages)
}
install_missing_packages(packages)
invisible(lapply(packages, library, character.only = TRUE))

## Simulation power estimate - Stop-level
#NOTE: We assume a 7.1% change in the arrest rates for stops related
#to the sex trades. We assume that we will have 800 treatment and 800
#comparison cases before the act and 1,000 treatment and 1,000 comparison 
#cases after the act. We explore sample sizes slightly higher and lower
#than that. And we assume a baseline arrest rate of 15.1%.

#Create a function that will simulate data with different sizes and baseline proportions
my_power_function <- function(baseline_prop, sample_sizetx1, sample_sizetx2, sample_sizec1, sample_sizec2) {
  #Store results here
  sig_results <- c()
  #Simulate data 1000 times
  for (i in 1:1000) {
    #Create the pre treatment group
    treatment_pre <- data.frame(condition=rep(1, sample_sizetx1), time=rep(0, sample_sizetx1), recid=rbinom(n=sample_sizetx1, size=1, prob=baseline_prop))
    #Create the post treatment group
    treatment_post <- data.frame(condition=rep(1, sample_sizetx2), time=rep(1, sample_sizetx2), recid=rbinom(n=sample_sizetx2, size=1, prob=baseline_prop-0.071))
    #Create the pre comparison group
    comparison_pre <- data.frame(condition=rep(0, sample_sizec1), time=rep(0, sample_sizec1), recid=rbinom(n=sample_sizec1, size=1, prob=baseline_prop))
    #Create the post comparison group
    comparison_post <- data.frame(condition=rep(0, sample_sizec2), time=rep(1, sample_sizec2), recid=rbinom(n=sample_sizec2, size=1, prob=baseline_prop))
    #Combine all conditions and time periods
    sample <- rbind(treatment_pre, treatment_post, comparison_pre, comparison_post)
    #Run a rudimentary model
    mylogit <- glm(recid ~ condition + time + condition:time, data = sample, family = "binomial")
    summary(mylogit)
    #Store whether the impact of the legislation is significant at an 0.05 level
    sig_results[i]  <- tidy(mylogit)$p.value[4] <= .05
  }
  #Store mean of significance
  sig_results %>%
    mean() %>%
    return()
}

#Try different sample sizes for each condition
tx1try <- c(600, 700, 800)
tx2try <- c(800, 900, 1000)
c1try <- c(600, 700, 800)
c2try <- c(800, 900, 1000)

#Store the findings 
power_levels_stops <- c()

#Run
for (i in 1:3) {
  power_levels_stops[i] <- my_power_function(0.151, tx1try[i], tx2try[i], c1try[i], c2try[i])
}
power_levels_stops

# Where do we cross 80%?
power_results_stops <- tibble(txpresample = tx1try,
                        power = power_levels_stops)
power_results_stops

#Plot the results
ggplot(power_results_stops, 
       aes(x = txpresample, y = power)) +
  geom_line(color = 'red', linewidth = 1.5) + 
  #Add a horizontal line at 80%
  geom_hline(aes(yintercept = .8), linetype = 'dashed') + 
  #Make it look nice
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'TX Sample Size Pre Act', y = 'Power') +
  ggtitle("Power with Different Sample Sizes")

## Simulation power estimate - Neighborhood Level
#NOTE: We assume 512 treatment block groups and 4,608 comparison 
#block groups measured monthly. We assume a starting level
#of crime 62.23 crimes with a standard deviation of 8.37.
#We assume crime will increase by 1 per year. 
#We also vary this slightly. To match existing research, we use an OLS 
#model with outcomes as a rate; although for the full 
#study we will explore whether a Poisson-based model is a better fit and
#if invited for the full proposal, will attempt to model that in this power
#analysis. Also, rather than look at an annual average, we will explore
#power analyses that allow for yearly rates and account for clustering. 
#Limited reseach in this area makes these sorts of power analysis models
#difficult to use.

#Create a function that will simulate data with different sizes and changes
my_power_function_2 <- function(changetry, sample_sizetx1, sample_sizetx2, sample_sizec1, sample_sizec2) {
  #Store results here
  sig_results <- c()
  #Simulate data 1000 times
  for (i in 1:1000) {
    #Create the pre treatment group
    treatment_pre <- data.frame(condition=rep(1, sample_sizetx1), time=rep(0, sample_sizetx1), crime=rnorm(n=sample_sizetx1, mean=62.23, sd=8.37))
    #Create the post treatment group
    treatment_post <- data.frame(condition=rep(1, sample_sizetx2), time=rep(1, sample_sizetx2), crime=rnorm(n=sample_sizetx2, mean=62.23+changetry, sd=8.37))
    #Create the pre comparison group
    comparison_pre <- data.frame(condition=rep(0, sample_sizec1), time=rep(0, sample_sizec1), crime=rnorm(n=sample_sizec1, mean=62.23, sd=8.37))
    #Create the post comparison group
    comparison_post <- data.frame(condition=rep(0, sample_sizec2), time=rep(1, sample_sizec2), crime=rnorm(n=sample_sizec2, mean=62.23, sd=8.37))
    #Combine all conditions and time periods
    sample <- rbind(treatment_pre, treatment_post, comparison_pre, comparison_post)
    #Run a rudimentary model
    myols <- lm(crime ~ condition + time + condition:time, data = sample)
    summary(myols)
    #Store whether the impact of the legislation is significant at an 0.05 level
    sig_results[i]  <- tidy(myols)$p.value[4] <= .05
  }
  #Store mean of significance
  sig_results %>%
    mean() %>%
    return()
}

#Try different changes
changetry <- c(1.75, 2.0, 2.25, 2.5)

#Store the findings 
power_levels_neighs <- c()

#Run
for (i in 1:4) {
  power_levels_neighs[i] <- my_power_function_2(changetry[i], 512, 512, 4608, 4608)
}
power_levels_neighs

# Where do we cross 80%?
power_results_neigh <- tibble(changes = changetry,
                              power = power_levels_neighs)
power_results_neigh

#Plot the results
ggplot(power_results_neigh, 
       aes(x = changes, y = power)) +
  geom_line(color = 'red', linewidth = 1.5) + 
  #Add a horizontal line at 80%
  geom_hline(aes(yintercept = .8), linetype = 'dashed') + 
  #Make it look nice
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Changes in Crime', y = 'Power') +
  ggtitle("Power with Different Changes in Crime")


knitr::stitch('SB-357 Power Analysis.r')
