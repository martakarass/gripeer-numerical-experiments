
#' @author
#' Marta Karas <mkaras2@jhu.edu>
#'
#' @description
#' Script to produce table with numerical experiments summaries: 
#' - regularization parameters, 
#' - execution time.

## Clear workspace
rm(list = ls())

## Set working directory path
wd <- "/Users/martakaras/Dropbox/_PROJECTS/OLD/LogisticPEER/gripeer-numerical-experiments"

library(Rmisc)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(data.table)

COL_NAMES <- c(
  "scenario_idx",
  "scenario_param_name",
  "scenario_param_val",
  "base_matrix_name",
  "n",
  "p",
  "i",
  "lambQ_Gripeer",
  "lambR_Gripeer",
  "lambQ_Gripeer_c_ridge",
  "exectime_Gripeer",
  "exectime_Gripeer_c_ridge"
)

## Read data: Scenario 1 -------------------------------------------------------

## Results file path
results_scenario1_A.path <- file.path(wd, "results", "results_scenario_1a.txt")
results_scenario1_B.path <- file.path(wd, "results", "results_scenario_1b.txt")
results_scenario1_C.path <- file.path(wd, "results", "results_scenario_1c.txt")

## Results zip-compressed file path
results_scenario1_A.zip.path <- file.path(wd, "results", "results_scenario_1a.txt.zip")
results_scenario1_B.zip.path <- file.path(wd, "results", "results_scenario_1b.txt.zip")
results_scenario1_C.zip.path <- file.path(wd, "results", "results_scenario_1c.txt.zip")

## If unzipped version does not exist, unzip
if (!file.exists(results_scenario1_A.path)){
  unzip(results_scenario1_A.zip.path, exdir = file.path(wd, "results"))
  message(paste0(results_scenario1_A.zip.path, " unzipped."))
}
if (!file.exists(results_scenario1_B.path)){
  unzip(results_scenario1_B.zip.path, exdir = file.path(wd, "results"))
  message(paste0(results_scenario1_B.zip.path, " unzipped."))
}
if (!file.exists(results_scenario1_C.path)){
  unzip(results_scenario1_C.zip.path, exdir = file.path(wd, "results"))
  message(paste0(results_scenario1_C.zip.path, " unzipped."))
}

results_scenario1_A <- as.data.frame(fread(results_scenario1_A.path, stringsAsFactors = FALSE))
results_scenario1_B <- as.data.frame(fread(results_scenario1_B.path, stringsAsFactors = FALSE))
results_scenario1_C <- as.data.frame(fread(results_scenario1_C.path, stringsAsFactors = FALSE))
results_scenario1_A$p <- 66
results_scenario1_B$p <- 198
results_scenario1_C$p <- 528

results_scenario1 <- rbind(
  results_scenario1_A,
  results_scenario1_B,
  results_scenario1_C
)

results_scenario1 <- 
  results_scenario1 %>%
  mutate(
    scenario_idx = 1,
    scenario_param_name = "dissimilarity",
    scenario_param_val = SC_mat_perm_val
  ) %>%
  tidyr::separate(SC_mat, into = c("base_matrix_name", "foo"), sep = "\\.")

dat_sc1 <- 
  results_scenario1[, COL_NAMES] %>%
  distinct()

str(dat_sc1)
# 'data.frame':	16000 obs. of  12 variables:
# $ scenario_idx            : num  1 1 1 1 1 1 1 1 1 1 ...
# $ scenario_param_name     : chr  "diss" "diss" "diss" "diss" ...
# $ scenario_param_val      : num  0 0 0 0 0 0 0 0 0 0 ...
# $ base_matrix_name        : chr  "A1" "A1" "A1" "A1" ...
# $ n                       : int  100 100 100 100 100 100 100 100 100 100 ...
# $ p                       : num  66 66 66 66 66 66 66 66 66 66 ...
# $ i                       : int  1 2 3 4 5 6 7 8 9 10 ...
# $ lambQ_Gripeer           : num  120 240 105 153 195 ...
# $ lambR_Gripeer           : num  0.114 0.187 0.153 0.11 0.16 ...
# $ lambQ_Gripeer_c_ridge   : num  0.263 0.507 0.301 0.241 0.378 ...
# $ exectime_Gripeer        : num  12.23 7.26 10.03 12.16 8.81 ...
# $ exectime_Gripeer_c_ridge: num  1.172 0.998 1.072 1.235 1.033 ...  



## Read data: Scenario 2 -------------------------------------------------------

## Results file path
results.path <- file.path(wd, "results", "results_scenario_2.txt")

## Results zip-compressed file path
results.zip.path <- file.path(wd, "results", "results_scenario_2.txt.zip")

## If unzipped version does not exist, unzip
if (!file.exists(results.path)){
  unzip(results.zip.path, exdir = file.path(wd, "results"))
  message(paste0(results.zip.path, " unzipped."))
}

results_scenario2 <- as.data.frame(fread(results.path, stringsAsFactors = FALSE))
results_scenario2$SC_mat_Anumber <- sapply(results_scenario2$SC_mat, function(val) substring(val, 1, 2))
results_scenario2$SC_mat_Aletter <- sapply(results_scenario2$SC_mat, function(val) substring(val, 3, 3))

SC_mat_Aletter2num.df <- data.frame(
  SC_mat_Aletter = c("a", "b", "c", "d"),
  SC_mat_neg_num = c(1,3,7,10)
)
results_scenario2 <- 
  results_scenario2 %>%
  left_join(SC_mat_Aletter2num.df, by = "SC_mat_Aletter")

results_scenario2 <- 
  results_scenario2 %>%
  mutate(
    scenario_idx = 2,
    scenario_param_name = "signs switched ncol",
    scenario_param_val = SC_mat_neg_num,
    base_matrix_name = SC_mat_Anumber,
    p = 66
  ) 

dat_sc2 <- 
  results_scenario2[, COL_NAMES] %>%
  distinct()

str(dat_sc2)
# 'data.frame':	12800 obs. of  12 variables:
# $ scenario_idx            : num  2 2 2 2 2 2 2 2 2 2 ...
# $ scenario_param_name     : chr  "signs switched ncol" "signs switched ncol" "signs switched ncol" "signs switched ncol" ...
# $ scenario_param_val      : num  1 1 1 1 1 1 1 1 1 1 ...
# $ base_matrix_name        : chr  "A1" "A1" "A1" "A1" ...
# $ n                       : int  100 100 100 100 100 100 100 100 100 100 ...
# $ p                       : num  66 66 66 66 66 66 66 66 66 66 ...
# $ i                       : int  1 2 3 4 5 6 7 8 9 10 ...
# $ lambQ_Gripeer           : num  174 162 105 124 357 ...
# $ lambR_Gripeer           : num  0.1162 0.185 0.1533 0.0988 0.1841 ...
# $ lambQ_Gripeer_c_ridge   : num  0.273 0.471 0.301 0.23 0.434 ...
# $ exectime_Gripeer        : num  35.5 27.6 27.4 35 17.6 ...
# $ exectime_Gripeer_c_ridge: num  2.75 2.31 2.79 3.01 2.48 ...



## Read data: Scenario 3 -------------------------------------------------------

## Results file path
results.path <- file.path(wd, "results", "results_scenario_3.txt")

## Results zip-compressed file path
results.zip.path <- file.path(wd, "results", "results_scenario_3.txt.zip")

## If unzipped version does not exist, unzip
if (!file.exists(results.path)){
  unzip(results.zip.path, exdir = file.path(wd, "results"))
  message(paste0(results.zip.path, " unzipped."))
}

results_scenario3 <- as.data.frame(fread(results.path, stringsAsFactors = FALSE))

results_scenario3 <- 
  results_scenario3 %>%
  mutate(
    scenario_idx = 3,
    scenario_param_name = "density ratio",
    scenario_param_val = scMatrixDerivValVec/100,
    p = 66
  ) %>%
  tidyr::separate(iScMatrixBaseName, into = c("base_matrix_name", "foo"), sep = "\\.")

dat_sc3 <- 
  results_scenario3[, COL_NAMES] %>%
  distinct()

str(dat_sc3)
# 'data.frame':	3600 obs. of  12 variables:
# $ scenario_idx            : num  3 3 3 3 3 3 3 3 3 3 ...
# $ scenario_param_name     : chr  "density ratio" "density ratio" "density ratio" "density ratio" ...
# $ scenario_param_val      : num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ base_matrix_name        : chr  "A1" "A1" "A1" "A1" ...
# $ n                       : int  100 100 100 100 100 100 100 100 100 100 ...
# $ p                       : num  66 66 66 66 66 66 66 66 66 66 ...
# $ i                       : int  1 2 3 4 5 6 7 8 9 10 ...
# $ lambQ_Gripeer           : num  0.362 1.323 0.248 0.119 1.878 ...
# $ lambR_Gripeer           : num  0.0789 0.169 0.149 0.1522 0.0443 ...
# $ lambQ_Gripeer_c_ridge   : num  0.263 0.507 0.301 0.241 0.378 ...
# $ exectime_Gripeer        : num  15.2 11.6 12.2 13 15.5 ...
# $ exectime_Gripeer_c_ridge: num  3.6 3.05 3.59 3.74 3.28 ...


## Combine data from 3 simulation scenarios ------------------------------------

dat_sc_all <- rbind(
  dat_sc1,
  dat_sc2,
  dat_sc3
)

# dat_sc_all_agg <- 
#   dat_sc_all %>%
#   group_by(
#     scenario_idx,
#     scenario_param_name,
#     scenario_param_val,
#     base_matrix_name,
#     n,
#     p
#   ) %>%
#   summarise(
#     ## mean
#     lambQ_Gripeer_mean            = round(mean(lambQ_Gripeer), 1),
#     lambR_Gripeer_mean            = round(mean(lambR_Gripeer), 1),
#     lambQ_Gripeer_c_ridge_mean    = round(mean(lambQ_Gripeer_c_ridge), 1),
#     exectime_Gripeer_mean         = round(mean(exectime_Gripeer), 1),
#     exectime_Gripeer_c_ridge_mean = round(mean(exectime_Gripeer_c_ridge), 1),
#     ## sd
#     lambQ_Gripeer_sd              = round(sd(lambQ_Gripeer), 1),
#     lambR_Gripeer_sd              = round(sd(lambR_Gripeer), 1),
#     lambQ_Gripeer_c_ridge_sd      = round(sd(lambQ_Gripeer_c_ridge), 1),
#     exectime_Gripeer_sd           = round(sd(exectime_Gripeer), 1),
#     exectime_Gripeer_c_ridge_sd   = round(sd(exectime_Gripeer_c_ridge), 1)
#   ) %>%
#   as.data.frame()
# 
# 
# dat_sc_all_agg_displ <- 
#   dat_sc_all_agg %>%
#   mutate(
#     lambQ_Gripeer            = paste0(sprintf('%.1f',lambQ_Gripeer_mean), " (", sprintf('%.1f',lambQ_Gripeer_sd), ")"),
#     lambR_Gripeer            = paste0(sprintf('%.1f',lambR_Gripeer_mean), " (", sprintf('%.1f',lambR_Gripeer_sd), ")"),
#     lambQ_Gripeer_c_ridge    = paste0(sprintf('%.1f',lambQ_Gripeer_c_ridge_mean), " (", sprintf('%.1f',lambQ_Gripeer_c_ridge_sd), ")"),
#     exectime_Gripeer         = paste0(sprintf('%.1f',exectime_Gripeer_mean), " (", sprintf('%.1f', exectime_Gripeer_sd), ")"),
#     exectime_Gripeer_c_ridge = paste0(sprintf('%.1f',exectime_Gripeer_c_ridge_mean), " (", sprintf('%.1f', exectime_Gripeer_c_ridge_sd), ")")
#   ) %>%
#   select(
#     scenario_idx,
#     base_matrix_name,
#     # scenario_param_name,
#     scenario_param_val,
#     n,
#     p,
#     lambQ_Gripeer,
#     lambR_Gripeer,
#     lambQ_Gripeer_c_ridge,
#     exectime_Gripeer
#    # exectime_Gripeer_c_ridge
#   ) %>%
#   arrange(
#     scenario_idx, 
#     base_matrix_name,
#     # scenario_param_name,
#     scenario_param_val,
#     n,
#     p
#     )
# 
# dim(dat_sc_all_agg_displ)
# head(dat_sc_all_agg_displ)
# 
# stargazer::stargazer(dat_sc_all_agg_displ, summary=FALSE, rownames=TRUE, font.size = "tiny")

rd_fact <- 2

dat_sc_all_agg <- 
  dat_sc_all %>%
  group_by(
    scenario_idx,
    # scenario_param_name,
    # scenario_param_val,
    base_matrix_name,
    n,
    p
  ) %>%
  summarise(
    ## mean
    lambQ_Gripeer_mean            = round(mean(lambQ_Gripeer), 1),
    lambR_Gripeer_mean            = round(mean(lambR_Gripeer), 1),
    lambQ_Gripeer_c_ridge_mean    = round(mean(lambQ_Gripeer_c_ridge), 1),
    exectime_Gripeer_mean         = round(mean(exectime_Gripeer), 1),
    exectime_Gripeer_c_ridge_mean = round(mean(exectime_Gripeer_c_ridge), 1),
    ## median
    lambQ_Gripeer_median            = round(median(lambQ_Gripeer), rd_fact),
    lambR_Gripeer_median            = round(median(lambR_Gripeer), rd_fact),
    lambQ_Gripeer_c_ridge_median    = round(median(lambQ_Gripeer_c_ridge), rd_fact),
    exectime_Gripeer_median         = round(median(exectime_Gripeer), rd_fact),
    exectime_Gripeer_c_ridge_median = round(median(exectime_Gripeer_c_ridge), rd_fact),
    ## sd
    lambQ_Gripeer_sd              = round(sd(lambQ_Gripeer), 1),
    lambR_Gripeer_sd              = round(sd(lambR_Gripeer), 1),
    lambQ_Gripeer_c_ridge_sd      = round(sd(lambQ_Gripeer_c_ridge), 1),
    exectime_Gripeer_sd           = round(sd(exectime_Gripeer), 1),
    exectime_Gripeer_c_ridge_sd   = round(sd(exectime_Gripeer_c_ridge), 1)
  ) %>%
  as.data.frame()


dat_sc_all_agg_displ <- 
  dat_sc_all_agg %>%
  mutate(
    lambQ_Gripeer            = paste0(sprintf('%.1f',lambQ_Gripeer_median)),
    lambR_Gripeer            = paste0(sprintf('%.1f',lambR_Gripeer_median)),
    lambQ_Gripeer_c_ridge    = paste0(sprintf('%.1f',lambQ_Gripeer_c_ridge_median)),
    exectime_Gripeer         = paste0(sprintf('%.1f',exectime_Gripeer_median)),
    exectime_Gripeer_c_ridge = paste0(sprintf('%.1f',exectime_Gripeer_c_ridge_median))
  ) %>%
  select(
    scenario_idx,
    base_matrix_name,
    # scenario_param_name,
    # scenario_param_val,
    n,
    p,
    lambQ_Gripeer,
    lambR_Gripeer,
    lambQ_Gripeer_c_ridge,
    exectime_Gripeer
    # exectime_Gripeer_c_ridge
  ) %>%
  arrange(
    scenario_idx, 
    base_matrix_name,
    # scenario_param_name,
    # scenario_param_val,
    n,
    p
  )

dim(dat_sc_all_agg_displ)
head(dat_sc_all_agg_displ)

stargazer::stargazer(dat_sc_all_agg_displ, summary=FALSE, rownames=TRUE, font.size = "tiny")

