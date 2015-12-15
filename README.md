# Data and scripts for pepBayes manuscript
### Contents
##### *main.R*
The repository root directory contains an R script *main.R* that loads data and creates figures used in the paper. Other resources in the repository are referenced from *main.R*. R's current working directory must be set such that the three **_path* variables at the head of *main.R* point to the appropriate data, figures, and lib folders in the repository. Currently the file paths are correctly set to be compatible with Mac OS X and most other Unix-like file systems. Several package dependencies are listed at the top of “main.R”, and commands for package installation are listed alongside. We specially note that the *pepBayes* R package currently requires a compiler implementing OpenMP (such as GNU GCC) to install from source. The script conducts the following general operations:

* Load libraries
* Load data
* Run pepBayes on data
* Plot calls
* Conduct ROC analysis
* Summarize simulation (optionally, run simulation)
* Plot supplementary figures

##### */data*
The two data sets used in the paper are saved as *p_rv.RData* (RV144) and *p_v3.RData* (Vax003). All patient clinical metadata except for treatment status has been stripped. *sim_output.RData* contains the simulation output we generated for the paper.

##### */lib*
The script *run_pepbayes.R* executes the pepBayes model on the data sets, as done in the manuscript. *load\_hiv\_anno.R* loads annotations for the HIV Env sequence. *utils.R* contains several functions used in the analysis and plots for our figures. *simulation.R* executes the simulation used in the paper. The execution of “simulation.R” is commented out in *main.R*, as it takes several days to execute. The output we generated is saved in “./data/sim_output.RData”.

##### */figures*
This directory contains the paper figures as we generated them. Running *main.R* should regenerate these figures (and in fact will replace them if steps are not taken to prevent overwrite).

