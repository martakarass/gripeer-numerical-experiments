
#' @author 
#' Marta Karas <mkaras2@jhu.edu>
#' 
#' @description 
#' Script to produce manuscript Figure 7.

## Clear workspace
rm(list = ls())

## Set working directory path 
wd <- "/Users/martakaras/Dropbox/_PROJECTS/LogisticPEER/gripeer-numerical-experiments"


## -----------------------------------------------------------------------------

## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_7")

library(Rmisc)
library(dplyr)
library(ggplot2)
select <- dplyr::select
filter <- dplyr::filter
summarize <- dplyr::summarize
mutate <- dplyr::mutate

## Results file path
results.path <- file.path(wd, "results", "results_scenario_2.txt")

## Results zip-compressed file path
results.zip.path <- file.path(wd, "results", "results_scenario_2.txt.zip")

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

dat <- as.data.frame(fread(results.path, stringsAsFactors = FALSE))

str(dat)
# data.frame':	844800 obs. of  8 variables:
#  $ SC_mat               : chr  "A1a.txt" "A1a.txt" "A1a.txt" "A1a.txt" ...
#  $ SC_mat_perm_val      : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ n                    : int  100 100 100 100 100 100 100 100 100 100 ...
#  $ i                    : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ b_idx                : int  1 2 3 4 5 6 7 8 9 10 ...
#  $ b_true               : num  0.727 -0.936 -2.429 1.772 -3.068 ...
#  $ b_est_griPEER        : num  0.676 0.698 -0.101 0.978 -0.667 ...
#  $ b_est_griPEER_c_ridge: num  0.348 0.442 -0.374 0.75 -0.406 ...
sort(unique(dat$SC_mat))
# [1] "A1a.txt" "A1b.txt" "A1c.txt" "A1d.txt" "A2a.txt" "A2b.txt" "A2c.txt" "A2d.txt"
# [9] "A3a.txt" "A3b.txt" "A3c.txt" "A3d.txt" "A4a.txt" "A4b.txt" "A4c.txt" "A4d.txt"
sort(unique(dat$n))
# [1] 100

# Define matrices names
dat$SC_mat_Anumber <- sapply(dat$SC_mat, function(val) substring(val, 1, 2))
dat$SC_mat_Aletter <- sapply(dat$SC_mat, function(val) substring(val, 3, 3))

## Melt data
data.m <- melt(dat, id.vars = setdiff(names(dat) , c("b_est_griPEER", "b_est_griPEER_c_ridge")))
names(data.m)[which(names(data.m) == "variable")] <- "est_method"
names(data.m)[which(names(data.m) == "value")] <- "b_est"
## Check the expected dimensions correctness
all(dim(data.m) == c(1689600, 10))

## Further aggregate data
SC_mat_Aletter2num.df <- data.frame(
  SC_mat_Aletter = c("a", "b", "c", "d"),
  SC_mat_neg_num = c(1,3,7,10)
)
data.m.g <- 
  data.m %>%
  filter(SC_mat_perm_val == 0) %>% # INTENTIONAL
  mutate(SE = (b_true - b_est)^2) %>%
  group_by(SC_mat_Anumber, SC_mat_Aletter, SC_mat_perm_val, est_method, i) %>%
  summarize(cnt = n(),
            MSEr = sum(SE)/sum((b_true^2))) %>% 
  left_join(SC_mat_Aletter2num.df, by = "SC_mat_Aletter") %>% 
  as.data.frame()
## Check the expected dimensions correctness
all(dim(data.m.g) == c(3200, 8))
data.m.g$est_method <- factor(as.character(data.m.g$est_method),
                              levels = c("b_est_griPEER", "b_est_griPEER_c_ridge"),
                              labels = c("griPEER","logistic ridge"))
SC_mat_Anumber.levels <- c("A1",
                           "A3", # intentional reordering
                           "A2", # intentional reordering
                           "A4")
SC_mat_Anumber.labels <- SC_mat.labels

data.m.g$SC_mat <- factor(as.character(data.m.g$SC_mat_Anumber),
                          levels = SC_mat_Anumber.levels,
                          labels = SC_mat_Anumber.labels)

## Further summarize data
summarySE_groupvars <- c("SC_mat_neg_num","SC_mat","est_method")
data.m.g.se <- summarySE(data.m.g, measurevar = "MSEr", groupvars = summarySE_groupvars)
data.m.g.se$n_p <- "n = 100, p = 66"


## Plot data
plt.base_size <- 11
pd <- position_dodge(0.05)

plt <- 
  ggplot(data.m.g.se,
         aes(x = SC_mat_neg_num, y = MSEr, color = est_method, 
             group = est_method)) + 
  geom_line(position=pd, size = 0.5, aes(linetype = est_method)) +
  geom_errorbar(aes(ymin=MSEr-se, ymax=MSEr+se), width=1, position=pd) +
  scale_x_continuous(breaks=c(1,3,7,10),
                     labels=as.integer(c(1,3,7,10))) +
  facet_grid(n_p ~ SC_mat, scales = "free_x") + 
  theme_bw(base_size = plt.base_size)  + 
  labs(y = "MSEr mean +/- se bars", 
       x = TeX("$k$ number of columns/rows with negative entries in $A^{true}$"), 
       title = "", 
       color = "Estimation method",
       linetype = "Estimation method") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values=c("blue", "red")) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(legend.position = c(.15, 0.4), legend.justification = c("right", "top"), legend.box.background = element_rect()) + 
  theme(legend.background = element_rect(colour = 'black', size = 0.5)) + 
  theme(legend.key.size = unit(0.4, "cm")) + 
  theme(legend.key.width = unit(1, "cm")) + 
  theme(legend.title = element_text(size = plt.base_size-1)) #+ 
plot(plt)


plt.path <- file.path(plt_save_dir, "sim2B_mser_mean_and_errorbar-1.jpeg")
ggsave(plt.path, plot = plt, height = 2.8, width = 10, units = "in")


