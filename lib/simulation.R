# Author: gimholte
###############################################################################
# CLUSTER COMMAND TO EXECUTE SCRIPT
# sbatch -n 1 -c 4 -t 2-0 --output=out.txt --wrap="R --no-save --no-restore < simulation.R"

#### SET WORKING DIRECTORY TO PEPBAYES SUPPLEMENT FOLDER 
setwd()
lib_path = "lib/"
data_path = "data/"

library(GenomicRanges)
library(truncnorm)
library(pepStat)
library(pepBayes)
source(paste0(lib_path, "utils.R"))

checkFile = "checkpoint.RData"
tempFile = "tempCheckpoint.RData"

if (file.exists(checkFile)) {
    load(checkFile)
    cat("Resuming after iteration", iter, "\n")
    iter <- iter + 1     # Important to avoid infinite loop
} else {
    n_peptide = 500
    n_subject = 30
    nt = 4
    prior_list_tmp = list(m0 = 4, nu_beta0 = 1,
            m1 = 2, nu_beta1 = 4, a = 1.5, b = 1.5,
            s_alpha = 4, lambda_alpha = 4,
            s_eps = 40, lambda_eps = 4,
            nu = 4)
    
    begin_time = Sys.time()
    set.seed(19287)
    Nsims = 500
    iter = 1                   # The first iteration
    sim_out = vector(mode = 'list', length = Nsims)
}

if (iter <= Nsims) {         
    # In case we were interrupted after
    # the last iteration
    for (iter in iter:Nsims) {            # Note _iter_:Nsims
        cat("Iteration", iter, "\n")
        ## Do iteration iter
        sim = simulatePairedPeptideSet(n_peptide, n_subject,
                n_rep = 3, n_active_pep = floor(n_peptide / 10),
                n_pep_per_pos = 5, prior_list_tmp, 
                weight_method = "all_separate")
        
        out_mcmc = pepBayesEcm(sim$pSet, sim$position_data, control_id = "pre",
                n_iter = 30, n_isamp = 30, n_threads = nt)
        
        out_ecm = pepBayesMcmc(sim$pSet, sim$position_data, control_id = "pre",
                n_samples = 1000, n_thin = 20, n_burn = 2000,
                n_threads = nt)
        
        sim_out[[iter]] = list(ppb_ecm = ppb(out_ecm), ppb_mcmc = ppb(out_mcmc),
                true_gamma = sim$gamma)
        ## Save the results of the iteration.  
        ## (By first saving to a temporary file
        ## and then renaming it into the checkpointing file, 
        ## we guard against being interrupted while saving.)
        save.image(tempFile)
        file.rename(tempFile, checkFile)
    }
}

cat(Sys.time() - begin_time, "\n")
save(sim_out, file = paste0(data_path, "sim_output.RData"))

if (file.exists(checkFile)) file.remove(checkFile)
