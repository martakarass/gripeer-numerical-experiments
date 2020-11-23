
#' @author 
#' Marta Karas <mkaras2@jhu.edu>
#' 
#' @description 
#' Script to produce manuscript Figure 4.

## Clear workspace
rm(list = ls())

## Set working directory path 
wd <- "/Users/martakaras/Dropbox/_PROJECTS/LogisticPEER/gripeer-numerical-experiments"


## -----------------------------------------------------------------------------

## Define directory to structural connectivity matrix objects 
sc_mat_dir <- file.path("data", "sc_matrices_used_scenario_2")
## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_4")

library(data.table)
library(ggplot2)
library(mdpeer)
library(latex2exp)

## Define ggplot2 options
sc.base_size <- 7

neg.col.n.vec <- c(1,4,7,10)

sim2_sc.plt <- function(pref.tmp, pref.latex.tmp = "", plt_prefix){
  
  sc_mat_objects     <- list.files(sc_mat_dir, full.names = TRUE)
  sc_mat_objects.tmp <- sort(sc_mat_objects[grepl(pref.tmp, sc_mat_objects)])
  pref.tmp.NO   <- paste0(pref.tmp,".txt")
  sc_mat_objects.tmp.NO <- sort(sc_mat_objects[grepl(pref.tmp.NO, sc_mat_objects)])
  sc_mat_objects.tmp <- sort(setdiff(sc_mat_objects.tmp, sc_mat_objects.tmp.NO))
  
  for (i in 1:length(sc_mat_objects.tmp)){ 
    ## Generate plot
    A_true_path <- sc_mat_objects.tmp[i]
    A_true <- as.matrix(fread(A_true_path))
    neg.col.n <- neg.col.n.vec[i]
    A_true_title <- TeX(paste0(pref.latex.tmp, '$^{true}$ with ', neg.col.n, ' neg-valued column/row(s)'))
    plt <- vizu.mat(A_true, title = A_true_title, base_size = sc.base_size) 
    plt <- plt + guides(fill=FALSE)
    ## Save plot
    plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", i, ".jpeg"))
    ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in") 
  }
}

sim2_sc.plt("A1", "$A_1$", "sim2_A_true_A1")
sim2_sc.plt("A3", "$A_2$", "sim2_A_true_A2")
sim2_sc.plt("A2", "$A_3$", "sim2_A_true_A3")
sim2_sc.plt("A4", "$A_4$", "sim2_A_true_A4")
