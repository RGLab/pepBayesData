# Run pepBayes model as MCMC on RV144 and Vax003 data sets 
#
# This script assumes the data have been loaded already.
# 
# Author: gimholte
###############################################################################

# settings for MCMC iterations
mcmc_burn = 10000
mcmc_thin = 20
mcmc_samp = 2000

# strings identifying control samples / pre-treatment samples 
vax_ctl_id = "PRE"
rv_ctl_id = "Pre"

# number of parallel threads to use
nt = 8
mcmc_rv = pepBayesMcmc(p_rv, position_data = pep_hxb2,
        control_id = rv_ctl_id, n_samples = mcmc_samp, n_thin = mcmc_thin,
        n_burn = mcmc_burn, n_threads = nt)

mcmc_rv_upd = pepBayesMcmc(p_rv_upd, position_data = pep_hxb2,
        control_id = rv_ctl_id, n_samples = mcmc_samp, n_thin = mcmc_thin,
        n_burn = mcmc_burn, n_threads = nt)

mcmc_v3 = pepBayesMcmc(p_v3, position_data = pep_hxb2,
        control_id = vax_ctl_id, n_samples = mcmc_samp, n_thin = mcmc_thin,
        n_burn = mcmc_burn, n_threads = nt)
