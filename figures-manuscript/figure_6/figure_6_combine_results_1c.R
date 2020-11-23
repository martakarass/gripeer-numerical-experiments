
#' @author
#' Marta Karas <mkaras2@jhu.edu>
#'
#' @description
#' Script to produce manuscript Figure 6.

## Clear workspace
rm(list = ls())

## Set working directory path
wd <- "/Users/martakaras/Dropbox/_PROJECTS/OLD/LogisticPEER/gripeer-numerical-experiments"

library(dplyr)
library(data.table)

## -----------------------------------------------------------------------------

# Create mapping based on one existing results file
map_base_fpath <- file.path(wd, "results", "results_scenario_1a.txt")
map_base_df <- 
  fread(map_base_fpath) %>%
  as.data.frame() 
map_df <- 
  map_base_df %>%
  select(SC_mat, SC_mat_perm_val) %>%
  distinct() %>%
  group_by(SC_mat) %>%
  arrange(SC_mat_perm_val) %>%
  mutate(diffidx = dplyr::row_number()) 

# Pull fules from dir with cluster-generated results
dir_results_scenario_1c.path <- paste0(wd, "/results/dir_results_scenario_1c")
fnames_results_scenario_1c <- list.files(dir_results_scenario_1c.path, full.names = TRUE)

df_list <- lapply(fnames_results_scenario_1c, function(fpath_tmp){
  f_tmp <- 
    fread(fpath_tmp) %>% 
    as.data.frame() %>%
    tidyr::separate(SC_mat, into = c("SC_mat", "c2", "c3"), sep = "_") %>%
    tidyr::separate(c3, into = c("diffidx", "c4"), sep = "\\.") %>%
    select(-c(c2, c4, SC_mat_perm_val)) %>%
    mutate(diffidx = as.integer(diffidx)) %>%
    mutate(SC_mat = paste0(SC_mat, ".txt")) %>%
    left_join(map_df, by = c("SC_mat", "diffidx")) %>%
    select(-diffidx) %>%
    dplyr::arrange(SC_mat, SC_mat_perm_val, n, i)
  f_tmp <- f_tmp[, names(map_base_df)]
})

length(df_list)

df_combined <- do.call("rbind", df_list)

# Write to file
out_fpath <- file.path(wd, "results", "results_scenario_1c.txt")
fwrite(as.data.table(df_combined), out_fpath)

