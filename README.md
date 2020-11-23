griPEER - numerical experiments
================

This repository contains code scripts and data needed to reproduce results (figures) presented in the Numerical experiments section of manuscript: Damian Brzyski, Marta Karas, Beau M. Ances, Mario Dzemidzic, Joaquin Goni, Timothy W. Randolph, Jaroslaw Harezlak, "Connectivity-Informed Adaptive Regularization for Generalized Outcomes". 

## License 

License of MATLAB code of griPEER method implementation is located at `MATLAB/software/griPEER_license.txt`. 

## Software usage notes 

- MATLAB toolbox `penalized` must be installed before using our software and running simulations code. To install the toolbox, visit corresponding [article website](https://www.jstatsoft.org/article/view/v072i06) and download the file linked next to "penalized.zip: MATLAB source code including replication materials". Next, unzip the zipped directory. Next, in MATLAB, change current folder to the unzipped toolbox directory and type: `install_penalized`. 

## Directory guide

### `data` 

Contains text files with pre-computed matrices used in simulation scenario. 
- `sc_matrices_used_scenario_1a` - directory with text files with pre-computed structural connectivity matrices used in simulation scenario 1 (part a). 
- `sc_matrices_used_scenario_1b_p198` - directory with text files with pre-computed structural connectivity matrices used in simulation scenario 1 (part b). 
- `sc_matrices_used_scenario_1c_p528` - directory with text files with pre-computed structural connectivity matrices used in simulation scenario 1 (part c). 
- `sc_matrices_used_scenario_2` - directory with text files with pre-computed structural connectivity matrices used in simulation scenario 2.
- `sc_matrices_used_scenario_3` - directory with text files with pre-computed structural connectivity matrices used in simulation scenario 3.
- `Z_cov.txt` - real data-based covariance matrix. 

### `figures-manuscript`

Contains directories named `figure_2`, ..., `figure_9` with 

- `R` code to generate plot(s) building up a manuscript figure,
- `JPEG` files of generated plot(s) building up a manucript figure, 

for manuscripf figures named Figure 2, ..., Figure 9, respectively. 

### `MATLAB`

#### `MATLAB/simulation-scripts` 

- `run_power_fdr.m` - MATLAB code to run power/FDR simulation. The results of this file execution are saved at `results/results_power_fdr.txt`.  
- `run_scenario_1a_p66.m` - MATLAB code to run simulation scenario 1 (part a). The results of this file execution are saved at `results/results_scenario_1a.txt`. 
- `run_scenario_1b_p198.m` - MATLAB code to run simulation scenario 1 (part b). The results of this file execution are saved at `results/results_scenario_1b.txt`.  
- `run_scenario_1c_p528_CLUSTER.m` - MATLAB code to run simulation scenario 1 (part c). The results of this file execution are saved at `results/results_scenario_1c.txt`.  
- `run_scenario_2.m` - MATLAB code to run simulation scenario 2. The results of this file execution are saved at `results/results_scenario_2.txt`.  
- `run_scenario_3.m` - MATLAB code to run simulation scenario 3. The results of this file execution are saved at `results/results_scenario_3.txt`.  
- `utils` - directory with auxiliary functions used in code to run simulations.

#### `MATLAB/software` 

- `griPEER.m` -  MATLAB code of griPEER method implementation. 
- `griPEER_license.txt` -  license of MATLAB code of griPEER method implementation. 

#### `MATLAB/software-usage-examples`

- `workingExample.m` - a working example of griPEER method code usage. 

### `results`

Contains text files with simulation results. 

- `results_power_fdr.txt` - text file with power/FDR simulation results. The results were generated with the code `MATLAB/simulation-scripts/run_power_fdr.m`.  
- `results_scenario_1a.txt` - text file with simulation scenario 1 (part a) results. The results were generated with the code `MATLAB/simulation-scripts/run_scenario_1a.m`. 
- `results_scenario_1b.txt` - text file with simulation scenario 1 (part b) results. The results were generated with the code `MATLAB/simulation-scripts/run_scenario_1b_p198.m`.  
- `results_scenario_1c.txt` - text file with simulation scenario 1 (part c) results. The results were generated with the code `MATLAB/simulation-scripts/run_scenario_1c_p528_CLUSTER.m` executed on a computing cluster. The bash files which run that  Maltab script are located at `bash_scripts/`.
- `results_scenario_2.txt` - text file with simulation scenario 2 results. The results were generated with the code `MATLAB/simulation-scripts/run_scenario_2.m`.   
- `results_scenario_3.txt` -text file with simulation scenario 3 results. The results were generated with the code `MATLAB/simulation-scripts/run_scenario_3.m`.   

### `bash_scripts`

Contains directories with bash files that run Maltab scripts for executing numerical experiment parts in a computing cluster. 

### `cluster_utils`

Contains utility files for executing numerical experiment parts in a computing cluster.


