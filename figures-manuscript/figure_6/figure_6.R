
#' @author
#' Marta Karas <mkaras2@jhu.edu>
#'
#' @description
#' Script to produce manuscript Figure 6.

## Clear workspace
rm(list = ls())

## Set working directory path
wd <- "/Users/martakaras/Dropbox/_PROJECTS/OLD/LogisticPEER/gripeer-numerical-experiments"


## -----------------------------------------------------------------------------

## Define directory to save plots
plt_save_dir <- file.path(wd, "figures-manuscript", "figure_6")

library(Rmisc)
library(dplyr)
library(ggplot2)
library(latex2exp)
select    <- dplyr::select
filter    <- dplyr::filter
summarize <- dplyr::summarize
mutate    <- dplyr::mutate

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


## -----------------------------------------------------------------------------

SC_mat.levels <- c("A1.txt",
                   "A3.txt", # intentional reordering
                   "A2.txt", # intentional reordering
                   "A4.txt")
SC_mat.labels <- c("base matrix: A1 \nhomologous regions",
                   "base matrix: A2 \nmodularity",
                   "base matrix: A3 \ndensity of connections, masked",
                   "base matrix: A4 \nneighboring regions")


## Read data results 1
results_A <- as.data.frame(fread(results_scenario1_A.path, stringsAsFactors = FALSE))

## Check expected aganist actual
## column names
col_names_exp <- c("SC_mat", "SC_mat_perm_val", "n", "i", "b_idx", "b_true", "b_est_griPEER", "b_est_griPEER_c_ridge")
col_names_actual <- names(results_A)[1:length(col_names_exp)]
all(col_names_exp == col_names_actual)
## structural connectivity matrices names
sc_names_exp <- c("A1.txt", "A2.txt", "A3.txt", "A4.txt")
sc_names_actual <- sort(unique(results_A$SC_mat))
all(sc_names_exp == sc_names_actual)
## nrow
nrow_exp <- 422400
nrow_actual <- nrow(results_A)
nrow_exp == nrow_actual
## n param
n_unique_exp <- c(100, 200)
n_unique_actual <- sort(unique(results_A$n))
all(n_unique_exp == n_unique_actual)


## Read data results 2
results_B <- as.data.frame(fread(results_scenario1_B.path, stringsAsFactors = FALSE))

## Check expected aganist actual
## column names
col_names_exp <- c("SC_mat", "SC_mat_perm_val", "n", "i", "b_idx", "b_true", "b_est_griPEER", "b_est_griPEER_c_ridge")
col_names_actual <- names(results_B)[1:length(col_names_exp)]
all(col_names_exp == col_names_actual)
## structural connectivity matrices names
sc_names_exp <- c("A1.txt", "A2.txt", "A3.txt", "A4.txt")
sc_names_actual <- sort(unique(results_B$SC_mat))
all(sc_names_exp == sc_names_actual)
## nrow
nrow_exp <- 633600
nrow_actual <- nrow(results_B)
nrow_exp == nrow_actual
## n param
n_unique_exp <- 100
n_unique_actual <- sort(unique(results_B$n))
n_unique_exp == n_unique_actual


## Read data results 3
results_C <- as.data.frame(fread(results_scenario1_C.path, stringsAsFactors = FALSE))

## Check expected aganist actual
## column names
col_names_exp <- c("SC_mat", "SC_mat_perm_val", "n", "i", "b_idx", "b_true", "b_est_griPEER", "b_est_griPEER_c_ridge")
col_names_actual <- names(results_C)[1:length(col_names_exp)]
all(col_names_exp == col_names_actual)
## structural connectivity matrices names
sc_names_exp <- c("A1.txt", "A2.txt", "A3.txt", "A4.txt")
sc_names_actual <- sort(unique(results_C$SC_mat))
all(sc_names_exp == sc_names_actual)
## nrow
nrow_exp <- 3379200
nrow_actual <- nrow(results_C)
nrow_exp == nrow_actual
## n param
n_unique_exp <- c(200, 400)
n_unique_actual <- sort(unique(results_C$n))
all(n_unique_exp == n_unique_actual)



# ------------------------------------------------------------------------------


## Add column with p value
results_A$p <- 66
results_B$p <- 198
results_C$p <- 528

## rbind two data frames. Reformat resulted data frame with results
dat <- rbind(results_A, results_B, results_C)
data.m <- reshape2::melt(dat, id.vars = setdiff(names(dat), c("b_est_griPEER", "b_est_griPEER_c_ridge")))
names(data.m)[which(names(data.m) == "variable")] <- "est_method"
names(data.m)[which(names(data.m) == "value")] <- "b_est"
data.m <- data.m %>% mutate(n_p = paste0(n, "_", p))
## Check dimmensions agreement with expected
# all(dim(data.m) == c(2112000, 10))
all(dim(data.m) == c(8870400, 16))

## Further mutate data
data.m.g <-
  data.m %>%
  mutate(SE = (b_true - b_est)^2) %>%
  group_by(SC_mat, SC_mat_perm_val, n, p, n_p, est_method, i) %>%
  summarize(cnt = n(),
            MSEr = sum(SE)/sum((b_true^2))) %>%
  as.data.frame()
## Check dimmensions agreement with expected
# all(dim(data.m.g) == c(19200, 9))
all(dim(data.m.g) == c(32000, 9))

# Change factor levels
sort(unique(data.m.g$n_p))
n_p_levels <- c(
  "100_66",
  "200_66",
  "100_198",
  "200_528",
  "400_528"
  )
n_p_labels <- c(
  "n = 100, p = 66",
  "n = 200, p = 66",
  "n = 100, p = 198",
  "n = 200, p = 528",
  "n = 400, p = 528"
  )
data.m.g$n_p <- factor(data.m.g$n_p,
                       levels = n_p_levels,
                       labels = n_p_labels)
data.m.g$est_method <- factor(as.character(data.m.g$est_method),
                              levels = c("b_est_griPEER", "b_est_griPEER_c_ridge"),
                              labels = c("griPEER","logistic ridge"))
data.m.g$SC_mat <- factor(as.character(data.m.g$SC_mat),
                          levels = SC_mat.levels,
                          labels = SC_mat.labels)

## Further summarize data (Gives count, mean, standard deviation, standard error of the
## mean, and confidence interval)
summarySE_groupvars <- c("SC_mat","SC_mat_perm_val","n_p","est_method")
data.m.g.se <- summarySE(data.m.g, measurevar = "MSEr", groupvars = summarySE_groupvars)

## Make and save the plot
pd <- position_dodge(0.05)
plt.base_size <- 11

plt <-
  ggplot(data.m.g.se,
         aes(x = SC_mat_perm_val, y = MSEr, color = est_method,
             group = est_method)) +
  geom_line(position=pd, size = 0.5, aes(linetype = est_method)) +
  geom_errorbar(aes(ymin=MSEr-se, ymax=MSEr+se), width=.05, position=pd) +
  facet_grid(n_p ~ SC_mat, scales = "free_x") +
  theme_bw(base_size = plt.base_size)  +
  labs(y = "MSEr mean +/- se bars",
       x = TeX("dissimilarity between $A^{obs}$ and $A^{true}$"),
       title = "",
       color = "Estimation method",
       linetype = "Estimation method") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values=c("blue", "red")) +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = c(.15, .123), legend.justification = c("right", "top"),
        legend.box.background = element_rect()) +
  theme(legend.background = element_rect(colour = 'black', size = 0.5)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  theme(legend.key.width = unit(1, "cm")) +
  theme(legend.title = element_text(size = plt.base_size-1))
plot(plt)

plt.path <- file.path(plt_save_dir, "sim1_mser_mean_and_errorbar-1.jpeg")
ggsave(plt.path, plot = plt, height = 9, width = 10, units = "in")


## Short plot version ----------------------------------------------------------

plt2 <-
  ggplot(data.m.g.se %>% filter(n_p == "n = 100, p = 66"),
         aes(x = SC_mat_perm_val, y = MSEr, color = est_method,
             group = est_method)) +
  geom_line(position=pd, size = 0.5, aes(linetype = est_method)) +
  geom_errorbar(aes(ymin=MSEr-se, ymax=MSEr+se), width=.05, position=pd) +
  facet_grid(n_p ~ SC_mat , scales = "free_x") +
  theme_bw(base_size = plt.base_size)  +
  labs(y = "MSEr mean +/- se bars",
       x = TeX("dissimilarity between $A^{obs}$ and $A^{true}$"),
       title = "",
       color = "Estimation method",
       linetype = "Estimation method") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_color_manual(values=c("blue", "red")) +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = c(.15, 0.4), legend.justification = c("right", "top"), legend.box.background = element_rect()) +
#theme(legend.position = c(.15, .123),
#       legend.justification = c("right", "top"),
#       legend.box.background = element_rect()) +
  theme(legend.background = element_rect(colour = 'black', size = 0.5)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  theme(legend.key.width = unit(1, "cm")) +
  theme(legend.title = element_text(size = plt.base_size-1))
plot(plt2)

plt2.path <- file.path(plt_save_dir, "sim1_mser_mean_and_errorbar_SHORT-1.jpeg")
ggsave(plt2.path, plot = plt2, height = 2.8, width = 10, units = "in")
