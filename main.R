# Script for generating pepBayes figures
# Author: gimholte
#############################################################

#### SET WORKING DIRECTORY TO PEPBAYES SUPPLEMENT FOLDER 
#### NOTE THAT FILE ADDRESS FOR *_path WILL DEPEND ON OPERATING SYSTEM
#setwd()
data_path = "data/"
fig_path = "figures/"
lib_path = "lib/"
#############################################################

#### LOAD LIBRARIES REQUIRED FOR ANALYSIS AND PLOTTING                 
#### COMMENTS INDICATE POSSIBLE INSTALLATION METHOD AND PACKAGE SOURCE
#### PACKAGE DEPENDENCIES ALSO LISTED (to our best knowledge)

# library(devtools)
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("biovizBase", "Biostrings", "GenomicRanges", "IRanges", "Gviz"))
# install.packages("data.table")
# install_github("RGLab/PViz")
# biocLite("pepStat")
# biocLite("pepDat")
# install.packages("Rcpp")
# install.packages("RcppArmadillo")
# install_github("RGLab/pepBayes")
# install.packages("ggplot2")
# install.packages("reshape2")
# install_github("RGLab/HIV.db")

library(Pviz)
library(pepStat)
library(pepDat)
data(pep_hxb2)
library(pepBayes)
library(ggplot2)
library(reshape2)
library(HIV.db)

#### LOAD DATA SETS
load(paste0(data_path, "p_rv.RData"))
load(paste0(data_path, "p_v3.RData"))

#### CREATE DUMMY UNPAIRED DATA SETS
p_rv_upd = p_rv
pData(p_rv_upd)$ptid = 1:200
p_v3_upd = p_v3
pData(p_v3_upd)$ptid = 1:180

#### SOURCE SCRIPT "run_pepbayes" WILL FIT MCMC PEPBAYES MODELS.

# RUNNING 3 PEPBAYES FITS WILL TAKE A WHILE (e.g. a few hours on a 2013 MacBook Pro), BUT
# THE SOURCED SCRIPT CONTAINS RUN PARAMETERS AND WILL FIT THE MODELS.
source(paste0(lib_path, "run_pepbayes.R"))
source(paste0(lib_path, "load_hiv_anno.R"))
source(paste0(lib_path, "utils.R"))

#################################################################
#### PLOT PROPORTION OF CALLS PER PEPTIDE
#################################################################

#### RV144 MCMC PAIRED, PLOT CALLS
rv144_trt_ind = mcmc_rv@slideInfo$metadata$treatment[seq(1, 200, 2)]
rv144_pos = mcmc_rv@slideInfo$orig_pos
freq_rv_mcmc = getFrequencyMatrix(mcmc_rv, rv144_trt_ind, fdr = .05)
limits = c(0, 1)
data_track_placebo = DTrack(start = rv144_pos, width = 1,
        data = freq_rv_mcmc$freq[1, , drop = FALSE],
        name = "Placebo, FDR = .05",
        legend = TRUE, alpha = .7, ylim = limits)

data_track_trt = DTrack(start = rv144_pos, width = 1,
        data = freq_rv_mcmc$freq[2, , drop = FALSE],
        name = "Treatment, FDR = .05",
        legend = TRUE, alpha = .7, ylim = limits)

pdf(paste0(fig_path, "calls_rv_mcmc.pdf"), width = 8, height = 5)
plotTracks(list(gp_track, loop_track[1:5], data_track_trt,
                data_track_placebo, axis_track), showFeatureID = TRUE,
        main = "RV144 response frequency with paired MCMC", col.axis = "black",
        fontcolor.title = "black", cex.main = 1)
dev.off()

#### RV144 MCMC UNPAIRED, PLOT CALLS
rv144_trt_ind_upd = with(mcmc_rv_upd@slideInfo$metadata, treatment[visit == "Post"])
rv144_pos_upd = mcmc_rv_upd@slideInfo$orig_pos
freq_rv_mcmc_upd = getFrequencyMatrix(mcmc_rv_upd,
        rv144_trt_ind_upd, fdr = .05)
limits = c(0, 1)
data_track_placebo_upd = DTrack(start = rv144_pos_upd, width = 1,
        data = freq_rv_mcmc_upd$freq[1, , drop = FALSE],
        name = "Placebo, FDR = .05",
        legend = TRUE, alpha = .7, ylim = limits)

data_track_trt_upd = DTrack(start = rv144_pos_upd, width = 1,
        data = freq_rv_mcmc_upd$freq[2, , drop = FALSE],
        name = "Treatment, FDR = .05",
        legend = TRUE, alpha = .7, ylim = limits)

pdf(paste0(fig_path, "calls_rv_mcmc_upd.pdf"), width = 8, height = 5)
plotTracks(list(gp_track, loop_track[1:5], data_track_trt_upd,
                data_track_placebo_upd, axis_track), showFeatureID = TRUE,
        main = "RV144 response frequency with unpaired MCMC", col.axis = "black",
        fontcolor.title = "black", cex.main = 1)
dev.off()


#### VAX003 MCMC PAIRED, PLOT CALLS
vax3_trt_ind = mcmc_v3@slideInfo$metadata$treatment[seq(1, 180, 2)]
freq_vax3_mcmc = getFrequencyMatrix(mcmc_v3, vax3_trt_ind, fdr =.05)
vax3_pos = mcmc_v3@slideInfo$orig_pos
limits = c(0, 1)

data_track_mcmc_vax3 = DTrack(start = vax3_pos[!(vax3_pos %in% c(548, 551))], width = 1,
        data = freq_vax3_mcmc$freq[, !(vax3_pos %in% c(548, 551))],
        name = "MCMC, FDR = .05",
        legend = TRUE, alpha = .7, ylim = limits)
pdf(paste0(fig_path, "calls_v3_mcmc.pdf"), width = 6, height = 3)
plotTracks(list(gp_track, loop_track[1:5], data_track_mcmc_vax3,
                axis_track), showFeatureID = TRUE,
        main = "Vax003 response frequency, paired MCMC", col.axis = "black",
        fontcolor.title = "black", cex.main = 1)
dev.off()


#################################################################
#### ROC ANALYSIS FOR PEPBAYES AND PEPSTAT METHODS
#################################################################

# peptide positions in gp41 used for both data sets
gp41_pos = c(512:900)

#### RV144 ROC
# difference in average fluorescence before and after treatmet
rv_postpre_diff = with(mcmc_rv@slideInfo, {
            # post-treatment minus pre-treatment fluorescence
            y_mean[, metadata$visit == "Post"] -
                y_mean[, metadata$visit == "Pre"]
        })

rv_truth = makeRocStandard(rv_postpre_diff, .01, gp41_pos,
        mcmc_rv@slideInfo$orig_pos,
        with(mcmc_rv@slideInfo$metadata, treatment[visit == "Post"]))

rv_truth_upd = makeRocStandard(rv_postpre_diff, .01, gp41_pos,
        mcmc_rv_upd@slideInfo$orig_pos,
        with(mcmc_rv_upd@slideInfo$metada, treatment[visit == "Post"]))
        
ROC_rv_pbayes = makeROC(ppb(mcmc_rv), rv_truth)
ROC_rv_pbayes_upd = makeROC(ppb(mcmc_rv_upd), rv_truth_upd)


#### VAX003 ROC
v3_postpre_diff = with(mcmc_v3@slideInfo, {
            y_mean[, metadata$visit == "POST"] -
                    y_mean[, metadata$visit == "PRE"]
        })

v3_truth = makeRocStandard(v3_postpre_diff, .01, gp41_pos,
        mcmc_v3@slideInfo$orig_pos,
        with(mcmc_v3@slideInfo$metadata, treatment[visit == "POST"]))

ROC_v3_pbayes = makeROC(ppb(mcmc_v3), v3_truth)
ROC_v3_pbayes_upd = makeROC(ppb(mcmc_v3_upd), v3_truth)

##### pepStat ROC
ps_v3 = summarizePeptides(p_v3, "median", pep_hxb2)
ps_rv = summarizePeptides(p_rv, "median", pep_hxb2)
pn_v3 = normalizeArray(ps_v3)
pn_rv = normalizeArray(ps_rv)

v3_post_idx = pData(pn_v3)$visit == "POST"
rv_post_idx = pData(pn_rv)$visit == "Post"

ROC_rv_pstat = makeROC(exprs(pn_rv)[,rv_post_idx] - exprs(pn_rv)[,!rv_post_idx],
        rv_truth)
ROC_v3_pstat = makeROC(exprs(pn_v3)[,v3_post_idx] - exprs(pn_v3)[,!v3_post_idx],
        v3_truth)


##### BIND ROCs INTO DATA FRAME FOR PLOTTING
ROC_df = rbind(ROC_rv_pbayes, ROC_rv_pbayes_upd, ROC_rv_pstat,
        ROC_v3_pbayes, ROC_v3_pbayes_upd, ROC_v3_pstat)
ROC_df$Dataset = rep(c("RV144", "Vax003"), each = 3 * nrow(ROC_rv_pbayes))
ROC_df$Method = rep(c("pepBayes paired", "pepBayes unpaired", "pepStat", 
                "pepBayes paired", "pepBayes unpaired", "pepStat"), 
        each = nrow(ROC_rv_pbayes))

pdf(paste0(fig_path, "roc_all.pdf"), width = 7.5, height = 3.5)
p = ggplot(ROC_df, aes(x = FPR, y = TPR)) + facet_grid(~Dataset) +
        geom_line(aes(linetype = Method)) + theme_bw() +
        xlab("False Positive Rate") +
        ylab("True Positive Rate") +
        ggtitle("Receiver Operating Characteristics")
print(p)
dev.off()

#################################################################
### SUMMARIZE SIMULATION RESULTS
#################################################################

# COMMENTED SOURCE FILE RUNS SIMULATION, WHICH WOULD TAKE SEVERAL DAYS
# ON A SINGLE DESKTOP COMPUTER. OUTPUT IS INCLUDED IN DATA FILES.
# source(paste0(lib_path, "simulation.R"))
load(paste0(data_path, "sim_output.RData"))
nreplicate = length(sim_out)

ecm_rates = array(0, dim = c(4, 2, nreplicate))
mcmc_rates = array(0, dim = c(4, 2, nreplicate))

for(i in 1:nreplicate) {
    ecm_rates[, , i] = simFDR(sim_out[[i]]$ppb_ecm, sim_out[[i]]$true_gamma)
    mcmc_rates[, , i] = simFDR(sim_out[[i]]$ppb_mcmc, sim_out[[i]]$true_gamma)
}

ecm_fdr = apply(ecm_rates, MARGIN = c(1, 2), mean)
mcmc_fdr = apply(mcmc_rates, MARGIN = c(1, 2), mean)
ecm_fdr_sd = apply(ecm_rates, MARGIN = c(1, 2), sd)
mcmc_fdr_sd = apply(mcmc_rates, MARGIN = c(1, 2), sd)

#################################################################
### SUPPLEMENTAL FIGURES IN WEB APPENDIX
#################################################################

#### WEB APPENDIX C
#### EXAMINE SLIDE EFFECTS IN PLACEBO DATA 
p_rv_sum = summarizePeptides(p_rv, summary = "mean", pep_hxb2)
pre_ind = pData(p_rv_sum)$visit == "Pre"
post_ind = pData(p_rv_sum)$visit == "Post"
plac_ind = pData(p_rv_sum)$treatment == "PLACEBO"

pre_dat = data.frame(exprs(p_rv_sum)[, pre_ind & plac_ind])
post_dat = data.frame(exprs(p_rv_sum)[, post_ind & plac_ind])

# SHOULD EVALUATE TO TRUE, INDICATING COLUMN ALIGNMENT
all(pData(p_rv_sum)$ptid[pre_ind] == pData(p_rv_sum)$ptid[post_ind])

# MATRICIES OF EXPRESSION DIFFERENCE AND EXPRESSION AVERAGE
diff_dat = melt(post_dat - pre_dat, variable.name = "Slide", 
        value.name = "Fluorescence_diff")
avg_dat = melt(.5 * (pre_dat + post_dat), variable.name = "Slide",
        value.name = "Fluorescence_avg")
avg_dat$Fluorescence_diff = diff_dat$Fluorescence_diff
levels(avg_dat$Slide) = paste0("Subject ", as.character(1:20))

# PLOT DIFFERENCE AGAINST AVERAGE EXPRESSION BY SUBJECT FOR EACH
# PLACEBO SUBJECT. LOESS SMOOTHER INDICATES TREND.
pdf(paste0(fig_path, "WA_C_additivity_loess.pdf"), width = 7, height = 8)
p = ggplot(avg_dat, aes(x = Fluorescence_avg, y = Fluorescence_diff)) +
        facet_wrap(~Slide, ncol = 4) + geom_point(alpha = .05) + 
        stat_smooth(method = "loess", enp.target = 12) +
        geom_hline(aes(yintercept = 0), colour = "red", lty = 2) +
        ggtitle("Peptide cross-slide expression differences\nPlacebo subjects, RV144 data")+
        xlab("Average peptide fluorescence, between slide pairs") +
        ylab("Difference of peptide fluorescence, between slide pairs") +
        theme_bw()
print(p)
dev.off()

#### WEB APPENDIX D
#### CREATE ROC CURVES AT ADDITIONAL ALPHA VALUES 

rv_truth_list = lapply(c(.01, .05, .1),
        makeRocStandard, reference_diff = rv_postpre_diff,
        neg_pos = gp41_pos,
        orig_pos = mcmc_rv@slideInfo$orig_pos,
        trt_factor = with(mcmc_rv@slideInfo$metadata, treatment[visit == "Post"]))

ROC_rv_list = lapply(rv_truth_list, makeROC, criteria = ppb(mcmc_rv))
ROC_rv_pstat_list = lapply(rv_truth_list, makeROC, criteria = rv_postpre_diff)

pdf(paste0(fig_path, "WA_D_roc_alpha.pdf"), width = 7, height = 4)
par(mfrow = c(1, 2))
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR",
        ylab = "TPR", main = "paired pepBayes ROC")
mapply(FUN = lines, ROC_rv_list, col = c("red", "blue", "black"))
legend("bottomright", bty = "n",
        leg = c(expression(textstyle(FWER) ~ alpha ~ "= .01"), 
                expression(textstyle(FWER) ~ alpha ~ "= .05"), 
                expression(textstyle(FWER) ~ alpha ~ "= .1")),
        lty = 1, col = c("red", "blue", "black"), lwd = 2)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "FPR",
        ylab = "TPR", main = "pepStat ROC")
mapply(FUN = lines, ROC_rv_pstat_list, col = c("red", "blue", "black"))
legend("bottomright", bty = "n",
        leg = c(expression(textstyle(FWER) ~ alpha ~ "= .01"), 
                expression(textstyle(FWER) ~ alpha ~ "= .05"), 
                expression(textstyle(FWER) ~ alpha ~ "= .1")),
        lty = 1, col = c("red", "blue", "black"), lwd = 2)
dev.off()


