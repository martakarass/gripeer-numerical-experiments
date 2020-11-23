
#' @author 
#' Marta Karas <mkaras2@jhu.edu>
#' 
#' @description 
#' Script to produce manuscript Figure 2. 

## Clear workspace
rm(list = ls())

## Set working directory path 
wd <- "/Users/martakaras/Dropbox/_PROJECTS/LogisticPEER/gripeer-numerical-experiments"


## -----------------------------------------------------------------------------

## Define directory to structural connectivity matrix objects 
sc_mat_dir <- file.path(wd, "data", "sc_matrices_used_scenario_1a")
## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_2")
## Prefix for plot output names
plt_prefix <- "sim_plt_A_true_base"

library(data.table)
library(ggplot2)
library(mdpeer)
library(latex2exp)

## Define ggplot2 options
sc.base_size <- 7

## Function to compute density of a matrix
dens.A <- function(A){
  sum(A != 0)/(dim(A)[1]*dim(A)[2])
} 

## Matrix objects path
pref.tmp   <- "A1_"
sc_mat.file.names <- paste0("A", 1:4, "_diffidx_1.txt")
sc_mat.file.paths <- file.path(sc_mat_dir, sc_mat.file.names)

## Matrix objects plot tile
sc_mat.labels <- c("$A_1$", "$A_3$", "$A_2$", "$A_4$")

for (i in 1:length(sc_mat.file.paths)){ 
  
  ## Generate plot
  sc_mat.file.path.i  <- sc_mat.file.paths[i]
  sc_mat.label.i      <- sc_mat.labels[i]
  A_obs               <- as.matrix(fread(sc_mat.file.path.i))
  A_obs               <- A_obs / max(A_obs)
  plt                 <- vizu.mat(A_obs, title = TeX(sc_mat.label.i), base_size = sc.base_size) 
  plt <- plt + guides(fill = FALSE)
  
  ## Save plot
  plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", i, ".jpeg"))
  ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in")
}
