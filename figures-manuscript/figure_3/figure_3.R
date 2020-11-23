
#' @author 
#' Marta Karas <mkaras2@jhu.edu>
#' 
#' @description 
#' Script to produce manuscript Figure 3.

## Clear workspace
rm(list = ls())

## Set working directory path 
wd <- "/Users/martakaras/Dropbox/_PROJECTS/LogisticPEER/gripeer-numerical-experiments"


## -----------------------------------------------------------------------------

## Define directory to structural connectivity matrix objects 
sc_mat_dir <- file.path(wd, "data", "sc_matrices_used_scenario_1a")
## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_3")

library(data.table)
library(ggplot2)
library(mdpeer)
library(latex2exp)

## Define ggplot2 options
sc.base_size <- 7


## -----------------------------------------------------------------------------
## Part 1 

pref.tmp   <- "A1_"
plt_prefix <- "sim1_A_obs_A1"

sc_mat_objects     <- list.files(sc_mat_dir, full.names = TRUE)
sc_mat_objects.tmp <- sort(sc_mat_objects[grepl(pref.tmp, sc_mat_objects)])
A_true_path        <- sc_mat_objects.tmp[1]
A_obs_paths.tmp    <- sc_mat_objects.tmp[-c(1)]
A_true             <- as.matrix(fread(A_true_path))
A_true             <- A_true / max(abs(A_true))
A_true.L           <- sum(abs(A_true)>0)

## Matrix "true"
plt_idx <- 1
A.true.plt.title  <- TeX('$A_1^{true}$: homologous regions') 
plt <- vizu.mat(A_true, title = A.true.plt.title, base_size = sc.base_size)
plt <- plt + guides(fill = FALSE) 
## Save plot
plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in")

## Matrix from scenario
for (A_obs_path in A_obs_paths.tmp){ # A_obs_path <- A_obs_paths.tmp[1]
  ## Increase "plot counter" used to generate file names
  plt_idx <- plt_idx + 1 
  ## Generate plot 
  A_obs <- as.matrix(fread(A_obs_path))
  A_obs <- A_obs / max(abs(A_obs))
  diss.denom.A <- abs(abs(A_true)-abs(A_obs))
  diss <- sum(diss.denom.A >0)/(A_true.L*2)
  plt.tile <- TeX(paste0("$A_1^{obs}$;  diss($A_1^{obs}$,$A_1^{true}$) = ", round(diss, 2)))
  plt <- vizu.mat(A_obs, title = plt.tile, base_size = sc.base_size) 
  plt <- plt + guides(fill = FALSE)
  ## Save plot
  if (plt_idx %in% c(1,3,6,8)){
    plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
    ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in") 
  }
}



## -----------------------------------------------------------------------------
## Part 2

pref.tmp   <- "A3_" ## "A3_" b/c  the naming of saved objects differs from order in manuscript
plt_prefix <- "sim1_A_obs_A2"

sc_mat_objects     <- list.files(sc_mat_dir, full.names = TRUE)
sc_mat_objects.tmp <- sort(sc_mat_objects[grepl(pref.tmp, sc_mat_objects)])
A_true_path        <- sc_mat_objects.tmp[1]
A_obs_paths.tmp    <- sc_mat_objects.tmp[-c(1)]
A_true             <- as.matrix(fread(A_true_path))
A_true             <- A_true / max(abs(A_true))
A_true.L           <- sum(abs(A_true)>0)

## Matrix "true"
plt_idx <- 1
A.true.plt.title  <- TeX('$A_2^{true}$: modularity') 
plt <- vizu.mat(A_true, title = A.true.plt.title, base_size = sc.base_size)
plt <- plt + guides(fill = FALSE) 
## Save plot
plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in")

## Matrix from scenario
for (A_obs_path in A_obs_paths.tmp){ # A_obs_path <- A_obs_paths.tmp[1]
  ## Increase "plot counter" used to generate file names
  plt_idx <- plt_idx + 1 
  ## Generate plot 
  A_obs <- as.matrix(fread(A_obs_path))
  A_obs <- A_obs / max(abs(A_obs))
  diss.denom.A <- abs(abs(A_true)-abs(A_obs))
  diss <- sum(diss.denom.A >0)/(A_true.L*2)
  plt.tile <- TeX(paste0("$A_2^{obs}$;  diss($A_2^{obs}$,$A_2^{true}$) = ", round(diss, 2)))
  plt <- vizu.mat(A_obs, title = plt.tile, base_size = sc.base_size) 
  plt <- plt + guides(fill = FALSE)
  ## Save plot
  if (plt_idx %in% c(1,3,6,8)){
    plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
    ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in") 
  }
}


## -----------------------------------------------------------------------------
## Part 3

pref.tmp   <- "A2_" ## "A2_" b/c  the naming of saved objects differs from order in manuscript
plt_prefix <- "sim1_A_obs_A3"

sc_mat_objects     <- list.files(sc_mat_dir, full.names = TRUE)
sc_mat_objects.tmp <- sort(sc_mat_objects[grepl(pref.tmp, sc_mat_objects)])
A_true_path        <- sc_mat_objects.tmp[1]
A_obs_paths.tmp    <- sc_mat_objects.tmp[-c(1)]
A_true             <- as.matrix(fread(A_true_path))
A_true             <- A_true / max(abs(A_true))
A_true.L           <- sum(abs(A_true)>0)

## Matrix "true"
plt_idx <- 1
A.true.plt.title  <- TeX('$A_3^{true}$: density of connections, masked') 
plt <- vizu.mat(A_true, title = A.true.plt.title, base_size = sc.base_size)
plt <- plt + guides(fill = FALSE) 
## Save plot
plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in")

## Matrix from scenario
for (A_obs_path in A_obs_paths.tmp){ # A_obs_path <- A_obs_paths.tmp[1]
  ## Increase "plot counter" used to generate file names
  plt_idx <- plt_idx + 1 
  ## Generate plot 
  A_obs <- as.matrix(fread(A_obs_path))
  A_obs <- A_obs / max(abs(A_obs))
  diss.denom.A <- abs(abs(A_true)-abs(A_obs))
  diss <- sum(diss.denom.A >0)/(A_true.L*2)
  plt.tile <- TeX(paste0("$A_3^{obs}$;  diss($A_3^{obs}$,$A_3^{true}$) = ", round(diss, 2)))
  plt <- vizu.mat(A_obs, title = plt.tile, base_size = sc.base_size) 
  plt <- plt + guides(fill = FALSE)
  ## Save plot
  if (plt_idx %in% c(1,3,6,8)){
    plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
    ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in") 
  }
}


## -----------------------------------------------------------------------------
## Part 4

pref.tmp   <- "A4_" 
plt_prefix <- "sim1_A_obs_A4"

sc_mat_objects     <- list.files(sc_mat_dir, full.names = TRUE)
sc_mat_objects.tmp <- sort(sc_mat_objects[grepl(pref.tmp, sc_mat_objects)])
A_true_path        <- sc_mat_objects.tmp[1]
A_obs_paths.tmp    <- sc_mat_objects.tmp[-c(1)]
A_true             <- as.matrix(fread(A_true_path))
A_true             <- A_true / max(abs(A_true))
A_true.L           <- sum(abs(A_true)>0)

## Matrix "true"
plt_idx <- 1
A.true.plt.title  <- TeX('$A_4^{true}$: neighboring regions') 
plt <- vizu.mat(A_true, title = A.true.plt.title, base_size = sc.base_size)
plt <- plt + guides(fill = FALSE) 
## Save plot
plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in")

## Matrix from scenario
for (A_obs_path in A_obs_paths.tmp){ # A_obs_path <- A_obs_paths.tmp[1]
  ## Increase "plot counter" used to generate file names
  plt_idx <- plt_idx + 1 
  ## Generate plot 
  A_obs <- as.matrix(fread(A_obs_path))
  A_obs <- A_obs / max(abs(A_obs))
  diss.denom.A <- abs(abs(A_true)-abs(A_obs))
  diss <- sum(diss.denom.A >0)/(A_true.L*2)
  plt.tile <- TeX(paste0("$A_4^{obs}$;  diss($A_4^{obs}$,$A_4^{true}$) = ", round(diss, 2)))
  plt <- vizu.mat(A_obs, title = plt.tile, base_size = sc.base_size) 
  plt <- plt + guides(fill = FALSE)
  ## Save plot
  if (plt_idx %in% c(1,3,6,8)){
    plt.path <- file.path(plt_save_dir, paste0(plt_prefix, "-", plt_idx, ".jpeg"))
    ggsave(plt.path, plot = plt, height = 2.4, width = 2.2, units = "in") 
  }
}


