
#' @author 
#' Marta Karas <mkaras2@jhu.edu>
#' 
#' @description 
#' Script to produce manuscript Figure 5.

## Clear workspace
rm(list = ls())

## Set working directory path 
wd <- "/Users/martakaras/Dropbox/_PROJECTS/LogisticPEER/gripeer-numerical-experiments"


## -----------------------------------------------------------------------------

## Define directory to structural connectivity matrix objects 
sc_mat_dir <- file.path("data", "sc_matrices_used_scenario_3")
## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_5")

library(data.table)
library(ggplot2)
library(mdpeer)
library(latex2exp)

sc.base_size <- 7
dens.new.to.orig.seq <- round(sort(unique(c(seq(50,100, length.out = 5), seq(100,150, length.out = 5)))))

plot.dens.rate <- function(pref.tmp, pref.latex.tmp = "", plt_prefix){
  A.names.tmp <- paste0(pref.tmp, "_", dens.new.to.orig.seq, ".txt")
  A.paths.tmp <- paste0(sc_mat_dir, "/", A.names.tmp)
  plt_idx <- 0
  for (A.path in A.paths.tmp){ 
    ## Generate plot
    plt_idx <- plt_idx + 1
    split.out <- strsplit(A.path, "_")[[1]]
    val <- split.out[length(split.out)]
    val <- round(as.numeric(strsplit(val, "\\.")[[1]][1])/100,2)
    A.tmp <- as.matrix(fread(A.path))
    A.tmp.title <- TeX(paste0(pref.latex.tmp, '$^{obs}$; density ratio = ', val))
    plt <- vizu.mat(A.tmp, title = A.tmp.title, base_size = sc.base_size) 
    plt <- plt + guides(fill=FALSE)
    ## Save the plot
    if (plt_idx %in% c(1,3,5,7,9)){
      plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
      ggsave(plt.path, plot = plt, height = 1.98, width = 1.88, units = "in") 
    }
  }
}

plot.dens.rate("A1", "$A_1$", "sim3_A_obs_A1")
plot.dens.rate("A3", "$A_2$", "sim3_A_obs_A2")
plot.dens.rate("A2", "$A_3$", "sim3_A_obs_A3")
plot.dens.rate("A4", "$A_4$", "sim3_A_obs_A4")


