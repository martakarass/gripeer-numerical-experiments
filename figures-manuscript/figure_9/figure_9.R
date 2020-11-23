
#' @author 
#' Marta Karas <mkaras2@jhu.edu>
#' 
#' @description 
#' Script to produce manuscript Figure 9.

## Clear workspace
rm(list = ls())

## Set working directory path 
wd <- "/Users/martakaras/Dropbox/_PROJECTS/LogisticPEER/gripeer-numerical-experiments"

## -----------------------------------------------------------------------------

## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_9")

library(Rmisc)
library(dplyr)
library(ggplot2)
library(data.table)
library(latex2exp)
select <- dplyr::select
filter <- dplyr::filter
summarize <- dplyr::summarize
mutate <- dplyr::mutate

## Results file path
results_power.path <- file.path(wd, "results", "results_power_fdr.txt") 

## Results zip-compressed file path
results_power.zip.path <- file.path(wd, "results", "results_power_fdr.txt.zip")

## If unzipped version does not exist, unzip
if (!file.exists(results_power.path)){
  unzip(results_power.zip.path, exdir = file.path(wd, "results"))
  message(paste0(results_power.zip.path, " unzipped."))
}



## -----------------------------------------------------------------------------

SC_mat.levels <- c("A1.txt",
                   "A3.txt", # intentional reordering
                   "A2.txt", # intentional reordering
                   "A4.txt")
SC_mat.labels <- c("base matrix: A1 \nhomologous regions",
                   "base matrix: A2 \nmodularity",
                   "base matrix: A3 \ndensity of connections, masked",
                   "base matrix: A4 \nneighboring regions")


## -----------------------------------------------------------------------------

## POWER

dat <- as.data.frame(fread(results_power.path, stringsAsFactors = FALSE))
dat$p <- 66
dat$n <- 150

dat.power <- dat %>%
  select(SC_mat, SC_mat_perm_val,n,i,POW_ASMP,POW_BOOT) %>%
  melt(id.vars = c("SC_mat", "SC_mat_perm_val",   "n", "i")) %>%
  tidyr::separate(col = variable, into = c("variable", "CI_type"), sep = "_") %>%
  mutate(n_p = paste0(n, "_", 66))

dat.fdr <- dat %>%
  select(SC_mat, SC_mat_perm_val,n,i,FDR_ASMP,FDR_BOOT) %>%
  melt(id.vars = c("SC_mat", "SC_mat_perm_val",   "n", "i")) %>%
  tidyr::separate(col = variable, into = c("variable", "CI_type"), sep = "_") %>%
  mutate(n_p = paste0(n, "_", 66))

dat.RelCoeff_frac <- dat %>%
  mutate(RelCoeff_frac = RelCoeff_n / 66) %>%
  select(SC_mat, SC_mat_perm_val,n,i,RelCoeff_frac) %>%
  melt(id.vars = c("SC_mat", "SC_mat_perm_val",   "n", "i")) %>%
  mutate(n_p = paste0(n, "_", 66))

dat.results <- rbind(dat.power, dat.fdr)
dat.results$CI_type <- factor(as.character(dat.results$CI_type),
                              levels = c("ASMP", "BOOT"),
                              labels = c("griPEER_asmp", "griPEER_boot"))

n.labels <- c("n = 150, p = 66")

## n_p
dat.results$n_p <- factor(dat.results$n_p, 
                          levels = sort(unique(dat.results$n_p)),
                          labels = n.labels)
dat.RelCoeff_frac$n_p <- factor(dat.RelCoeff_frac$n_p, 
                                levels = sort(unique(dat.RelCoeff_frac$n_p)),
                                labels = n.labels)

## SC_mat
dat.results$SC_mat <- factor(as.character(dat.results$SC_mat),
                             levels = SC_mat.levels,
                             labels = SC_mat.labels)
dat.RelCoeff_frac$SC_mat <- factor(as.character(dat.RelCoeff_frac$SC_mat),
                                   levels = SC_mat.levels,
                                   labels = SC_mat.labels)


dat.results.se <- Rmisc::summarySE(dat.results, 
                                   measurevar = "value", 
                                   groupvars=c("SC_mat","SC_mat_perm_val","n_p","variable", "CI_type"))
dat.RelCoeff_frac.se <- Rmisc::summarySE(dat.RelCoeff_frac, 
                                         measurevar = "value", 
                                         groupvars=c("SC_mat","SC_mat_perm_val","n_p","variable"))

plt.base_size <- 8
pd <- position_dodge(0.05)
plt.title <- ""

plt <- 
  ggplot(dat.results.se %>% filter(variable == "POW", SC_mat == unique(dat.results.se$SC_mat)[1]),  # unique(dat.results.se$SC_mat)[3]),
         aes(x = SC_mat_perm_val, 
             y = value, 
             color = CI_type)) + 
  geom_line(aes(linetype = CI_type), position=pd, size = 0.5) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.05, position=pd, linetype = 1,
                show.legend = FALSE) +
  facet_grid(n_p ~ SC_mat, scales = "free_x") + 
  theme_bw(base_size = plt.base_size)  + 
  labs(y = "power mean +/- se bars", 
       x = TeX("dissimilarity between $A^{obs}$ and $A^{true}$"), 
       title = plt.title, 
       color = "Method", 
       linetype = "Method") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values=c("green4", "red4")) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"), legend.box.background = element_rect()) + 
  theme(legend.background = element_rect(colour = 'black', size = 0.5)) + 
  theme(legend.key.size = unit(0.4, "cm")) + 
  theme(legend.key.width = unit(1, "cm")) + 
  theme(legend.title = element_text(size = plt.base_size-1))
plot(plt)

plt.path <- file.path(plt_save_dir, "sim_CI_mser_mean_and_errorbar_POWER_1mat-1.jpeg")
ggsave(plt.path, plot = plt, height = 3, width = 3, units = "in")



## -----------------------------------------------------------------------------

## FDR

plt <- 
  ggplot(dat.results.se %>% filter(variable == "FDR",  SC_mat == unique(dat.results.se$SC_mat)[1]),  # unique(dat.results.se$SC_mat)[3]),
         aes(x = SC_mat_perm_val, 
             y = value, 
             color = CI_type)) + 
  geom_line(aes(linetype = CI_type), position=pd, size = 0.5) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.05, position=pd, linetype = 1,
                show.legend = FALSE) +
  facet_grid(n_p ~ SC_mat, scales = "free_x") + 
  theme_bw(base_size = plt.base_size)  + 
  labs(y = "FDR mean +/- se bars", 
       x = TeX("dissimilarity between $A^{obs}$ and $A^{true}$"), 
       title = plt.title, 
       color = "Method",
       linetype = "Method") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values=c("green4", "red4")) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"), legend.box.background = element_rect()) + 
  theme(legend.background = element_rect(colour = 'black', size = 0.5)) + 
  theme(legend.key.size = unit(0.4, "cm")) + 
  theme(legend.key.width = unit(1, "cm")) + 
  theme(legend.title = element_text(size = plt.base_size-1)) 
plot(plt)


plt.path <- file.path(plt_save_dir, "sim_CI_mser_mean_and_errorbar_FDR_1mat-1.jpeg")
ggsave(plt.path, plot = plt, height = 3, width = 3, units = "in")

