
#' @author 
#' Marta Karas <mkaras2@jhu.edu>
#' 
#' @description 
#' Script to produce manuscript Figure 8.

## Clear workspace
rm(list = ls())

## Set working directory path 
wd <- "/Users/martakaras/Dropbox/_PROJECTS/OLD/LogisticPEER/gripeer-numerical-experiments"


## -----------------------------------------------------------------------------

## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_8")

library(Rmisc)
library(dplyr)
library(ggplot2)
select <- dplyr::select
filter <- dplyr::filter
summarize <- dplyr::summarize
mutate <- dplyr::mutate

## Results file path
results.path <- file.path(wd, "results", "results_scenario_3.txt")

## Results zip-compressed file path
results.zip.path <- file.path(wd, "results", "results_scenario_3.txt.zip")

## If unzipped version does not exist, unzip
if (!file.exists(results.path)){
  unzip(results.zip.path, exdir = file.path(wd, "results"))
  message(paste0(results.zip.path, " unzipped."))
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

## Read data with results, reshape
dat <- as.data.frame(fread(results.path, stringsAsFactors = FALSE))
data.m <- reshape2::melt(dat, id.vars = setdiff(names(dat) , c("b_est_griPEER", "b_est_griPEER_c_ridge")))
names(data.m)[which(names(data.m) == "variable")] <- "est_method"
names(data.m)[which(names(data.m) == "value")] <- "b_est"
## Check dimension correctness
# all(dim(data.m) == c(475200, 8))
all(dim(data.m) == c(475200, 14))

## Aggregate data
data.m.g <- 
  data.m %>%
  mutate(SE = (b_true - b_est)^2) %>%
  group_by(iScMatrixBaseName, scMatrixDerivValVec, est_method, i) %>%
  summarize(cnt = n(),
            MSEr = sum(SE)/sum((b_true^2))) %>% 
  as.data.frame()
all(dim(data.m.g) == c(7200, 6))

data.m.g$est_method <- factor(as.character(data.m.g$est_method),
                              levels = c("b_est_griPEER", "b_est_griPEER_c_ridge"),
                              labels = c("griPEER", "logistic ridge"))
data.m.g$iScMatrixBaseName <- factor(as.character(data.m.g$iScMatrixBaseName),
                                     levels = SC_mat.levels,
                                     labels = SC_mat.labels)

## Further summarize the data 
data.m.g.se <- Rmisc::summarySE(data.m.g, measurevar = "MSEr", groupvars=c("scMatrixDerivValVec","iScMatrixBaseName","est_method"))
data.m.g.se$n_p <- "n = 100, p = 66"

## Plot
plt.base_size <- 11
pd <- position_dodge(0.05)

plt <- 
  ggplot(data.m.g.se,
         aes(x = scMatrixDerivValVec/100, y = MSEr, color = est_method, 
             group = est_method)) + 
  geom_vline(xintercept = 1, linetype = 2, size = 1.5, color = "green", alpha = 0.6) + 
  geom_line(position=pd, size = 0.5, aes(linetype = est_method)) +
  geom_errorbar(aes(ymin=MSEr-se, ymax=MSEr+se), width=.05, position=pd) +
  facet_grid(n_p ~ iScMatrixBaseName, scales = "free_x") + 
  theme_bw(base_size = plt.base_size)  + 
  labs(y = "MSEr mean +/- se bars", 
       x = TeX("density ratio: $dens(A^{obs})/dens(A^{true})$ "), 
       title = "", 
       color = "Estimation method", 
       linetype = "Estimation method") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values=c("blue", "red")) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(legend.position = c(.15, .4), legend.justification = c("right", "top"), legend.box.background = element_rect()) + 
  theme(legend.background = element_rect(colour = 'black', size = 0.5)) + 
  theme(legend.key.size = unit(0.4, "cm")) + 
  theme(legend.key.width = unit(1, "cm")) + 
  theme(legend.title = element_text(size = plt.base_size-1))
plot(plt)

plt.path <- file.path(plt_save_dir, "sim3_mser_mean_and_errorbar-1.jpeg")
ggsave(plt.path, plot = plt, height = 2.8, width = 10, units = "in")


