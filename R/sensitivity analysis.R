TIME_CUT <- '2020-01-20'

source('R/utils.R')
source('R/facet_wrap_custom.R')
source("R/data preparation.R")
source("R/statistical_test_new.R")
test_results

TIME_CUT <- '2020-02-01'

source('R/utils.R')
source('R/facet_wrap_custom.R')
source("R/data preparation.R")
source("R/statistical_test_new.R")
test_results


library(future.apply)
sensitivity_parent <- function(seed){
  list2env(list(seed = seed, TIME_CUT = '2020-01-23'), 
           envir = environment())
  options(warn=-1)
  source('R/utils.R', local = TRUE)
  source("R/data preparation.R", local = TRUE)
  source("R/statistical_test_new.R", local = TRUE)
  return(test_results)
}

plan(multisession, workers = 8)

parent_results <- future_lapply(1:100,
                                FUN = sensitivity_parent,
                                future.packages = c("data.table", "dplyr", "network", "boot"),
                                future.seed = TRUE)

a <- rowSums(sapply(parent_results, function(x) {(x$p_value <= 0.05)})) / 100

parent_results[[1]]$measures[which(! (a<=0.05 | a>=0.95))]

results <- data.table(measures = parent_results[[1]]$measures,
                      significant_prop = a)
write.csv(results, "result/result_sample_parent.csv")
