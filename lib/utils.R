# Some functions used in creating pepBayes figures and simulation
# 
# Author: gimholte
###############################################################################

getFrequencyMatrix = function(pbayes_fit, group_indicator, fdr) {
    calls = pepBayesCalls(pbayes_fit, fdr = fdr)
    groups = unique(as.factor(group_indicator))
    ngroups = length(groups)
    freq_matrix = matrix(0, nrow = ngroups, ncol = nrow(calls))
    for(i in 1:length(groups)) {
        freq_matrix[i,] = rowMeans(calls[, group_indicator == groups[i]])        
    }
    return(list(freq = freq_matrix, groups = groups))
}

makeRocStandard = function(reference_diff, alpha, neg_pos, orig_pos,
        trt_factor) {
    m = rowMeans(reference_diff)
    n = ncol(reference_diff)
    p = nrow(reference_diff)
    t_rv = sqrt(n) *  m / apply(reference_diff, 1, sd)
    pval_rv = pt(t_rv, df = n - 1, lower.tail = FALSE)
    pos_idx = which(pval_rv < alpha / p)
    neg_idx = which(orig_pos %in% neg_pos)
    
    standard = matrix(FALSE, nrow = length(pos_idx) + length(neg_idx),
            ncol = n)
    vacc_idx = toupper(as.character(trt_factor)) == "VACCINE"
    standard[1:length(pos_idx), vacc_idx] = TRUE
    attr(standard, "pos_idx") = pos_idx
    attr(standard, "neg_idx") = neg_idx
    return(standard)
}

# not optimized but it works:
makeROC = function(criteria, standard) {
    sub_idx = c(attr(standard, "pos_idx"), attr(standard, "neg_idx"))
    criteria = criteria[sub_idx,]
    if (is.null(criteria) | is.null(standard))
        stop("non-null criteria and standard required")
    if (any(dim(criteria) != dim(standard)))
        stop("criteria and standard must have matching dimensions")
    quantile_vec = quantile(criteria, probs = seq(.999, .001, by = -.005),
            names = FALSE)
    ngrid = length(quantile_vec)
    ROC = sapply(1:ngrid, function (i) {
        calls = (criteria >= quantile_vec[i])
        TP = sum(calls & standard) 
        FN = sum(!calls & standard)
        TN = sum(!calls & !standard)
        FP = sum(calls & !standard)
        TPR = TP/(TP + FN)
        FPR = FP/(FP + TN)
        return(c(FPR, TPR))
    })
    ROC = t(cbind(c(0, 0), ROC))
    colnames(ROC) = c("FPR", "TPR")    
    return(as.data.frame(ROC))
}

createPeptideSequence = function(i) {
    base = 19
    seq = c()
    AA = c("G", "V", "L", "I", "M", "F", "W", "P",
            "S", "T", "C", "Y", "N", "Q", "D", "E", "K", "R", "H")
    while (i > 0) {
        bit = i %% base
        i = floor(i / base)
        seq = c(AA[bit + 1], seq)
    }
    
    len = length(seq)
    rest = rep("A", 15 - len)
    return(paste(c(rest, seq), collapse = ""))
}

createWeights = function(method, n_peptide, n_subject, nu) { 
    if (method == "normal") {
        w0 = matrix(1, n_peptide, n_subject)
        w1 = matrix(1, n_peptide, n_subject)
        u0 = matrix(1, n_peptide, n_subject)
        u1 = matrix(1, n_peptide, n_subject)
        return(list(w0 = w0, w1 = w1, u0 = u0, u1 = u1))
    }
    
    w0 = matrix(rgamma(n_subject * n_peptide, shape = nu / 2, 
                    rate = (nu - 2) / 2), nrow = n_peptide, ncol = n_subject)
    w1 = matrix(rgamma(n_subject * n_peptide, shape = nu / 2, 
                    rate = (nu - 2) / 2), nrow = n_peptide, ncol = n_subject)
    u0 = matrix(rgamma(n_subject * n_peptide, shape = nu / 2, 
                    rate = (nu - 2) / 2), nrow = n_peptide, ncol = n_subject)
    u1 = matrix(rgamma(n_subject * n_peptide, shape = nu / 2, 
                    rate = (nu - 2) / 2), nrow = n_peptide, ncol = n_subject)
    
    if (method == "shared") {
        u0 = w0
        u1 = w1
    }
    
    if (!any(method %in% c("shared", "all_separate", "normal"))) {
        stop("invalid weightings method")
    }
    return(list(w0 = w0, w1 = w1, u0 = u0, u1 = u1))
}

simFDR = function(ppa_out, standard, fdr_breaks = c(.01, .05, .1, .2)) {
    o = order(ppa_out, decreasing = TRUE)
    csum_order = cumsum(ppa_out[o])
    n_call = length(ppa_out)
    n_break = length(fdr_breaks)
    output = matrix(0, n_break, 2)
    fdr_vec = csum_order / (1:n_call)
    for (j in 1:n_break) {
        ind = which(fdr_vec > 1 - fdr_breaks[j])
        if (length(ind) == 0)
            calls = ppa_out > 1
        else
            calls = ppa_out > ppa_out[o][max(ind)]
        FP = sum((calls & (!standard)) & (standard == FALSE))
        TP = sum((calls & standard) & (standard == TRUE))
        total_positives = sum(standard)
        if (TP + FP == 0)
            output[j, 1] = 0
        else
            output[j, 1] = FP / (TP + FP)
        output[j, 2] = TP / total_positives       
    }
    return(output)
}

simulatePairedPeptideSet = function(n_peptide, n_subject,
        n_rep = 3, n_active_pep = floor(n_peptide / 10),
        n_pep_per_pos = 5, prior_list, weight_method) {
    n_inactive_pep = n_peptide - n_active_pep
    n_slide = 2 * n_subject
    
    gamma = matrix(0, nrow = n_peptide, ncol = n_subject)
    if (n_active_pep > 0) {
        active_pep_ind = rep(FALSE, n_peptide)
        active_pep_ind[sample(1:n_peptide, n_active_pep, replace = FALSE)] = TRUE
        active_prop = with(prior_list,
                rbeta(n_active_pep, shape1 = a, shape2 = b))
        beta1 = with(prior_list,
                {
                    b = rep(0, n_peptide)
                    b[active_pep_ind] =  rtruncnorm(n_active_pep,
                            a = 0, b = Inf, mean = m1, sd = 1 / sqrt(nu_beta1))
                    b
                })
        gamma[active_pep_ind,] = replicate(n_subject, 
                {
                    rbinom(n_active_pep, size = 1, prob = active_prop)
                })
    } else {
        beta1 = rep(0, n_peptide)
    }
    
    pre_ind = rep(c(TRUE, FALSE), n_subject)
    post_ind = !pre_ind
    
    beta0 = with(prior_list, 
            rnorm(n_peptide, mean = m0, sd = 1 / sqrt(nu_beta0)))
    
    mu = rnorm(n_slide, mean = 0, sd = 2)
    mu = mu - mean(mu)
    
    nu_eps = with(prior_list,
            rgamma(n_peptide, shape = s_eps, rate = lambda_eps))
    nu_alpha = with(prior_list,
            rgamma(n_peptide, shape = s_alpha, rate = lambda_alpha))    
    
    weights = createWeights(weight_method, n_peptide, n_subject, 
            prior_list$nu)
    
    w0 = weights$w0
    w1 = weights$w1
    u0 = weights$u0
    u1 = weights$u1
    
    alpha0_tmp = matrix(rnorm(n_subject * n_peptide,
                    mean = 0, sd = 1 / sqrt(u0 * nu_alpha)),
            nrow = n_peptide, ncol = n_subject)
    alpha0 = matrix(0, n_peptide, n_slide)
    alpha0[,pre_ind] = alpha0_tmp
    alpha0[,post_ind] = alpha0_tmp
    
    alpha1 = matrix(rnorm(n_subject * n_peptide,
                    mean = 0, sd = 1 / sqrt(u1 * nu_alpha)),
            nrow = n_peptide, ncol = n_subject)
    
    expr_mean = alpha0 + beta0 + rep(mu, each = n_peptide)
    expr_mean[,post_ind] = expr_mean[,post_ind] + 
            gamma * (beta1 + alpha1 - alpha0[,post_ind])
    
    expr = lapply(1:n_peptide, function(p) 
            {
                exp_mat = matrix(0, n_rep, n_slide)
                pre_weight = w0[p,]
                post_weight = w1[p,]
                if (weight_method == "shared") {
                    post_weight = gamma[p,] * w1[p,] +
                            (1 - gamma[p,]) * w0[p,]
                }
                exp_mat[,pre_ind] = t(replicate(n_rep,
                                rnorm(n_subject, mean = expr_mean[p, pre_ind],
                                        sd = 1 / sqrt(pre_weight * nu_eps[p]))))
                exp_mat[,post_ind] = t(replicate(n_rep,
                                rnorm(n_subject, mean = expr_mean[p, post_ind],
                                        sd = 1 / sqrt(post_weight * nu_eps[p]))))  
                return(exp_mat)
            })
    my_desc <- new("MIAME")
    expr = do.call(rbind, expr)
    slide_metadata = data.frame(ptid = rep(1:n_subject, each = 2),
            visit = rep(c("pre", "post"), n_subject))
    colnames(expr) = with(slide_metadata,
            paste(ptid, visit, sep = "_"))
    
    feature_sequence = rep(sapply(1:n_peptide, createPeptideSequence), 
            each = n_rep)
    feature_id = feature_sequence
    
    n_pos = ceiling(n_peptide / n_pep_per_pos)
    position = rep(1:n_pos, each = n_pep_per_pos)[1:n_peptide]
    position_data = data.frame(peptide = 
                    feature_sequence[seq(1, n_peptide * n_rep, by = n_rep)],
            start = position, width = 15, end = position + 14)
    
    pep_set = new('peptideSet',
            featureRange = GRanges(seqnames = " ", strand = "*",
                    ranges = IRanges(rep(0, n_peptide * n_rep),
                            rep(0, n_peptide * n_rep)),
                    featureID = feature_id,
                    peptide = feature_sequence),
            exprs = expr,
            experimentData = my_desc)
    rownames(slide_metadata) = colnames(expr)
    pData(pep_set) = slide_metadata
    
    return(list(pSet = pep_set,
                    position_data = position_data,
                    expr_mean = expr_mean, 
                    gamma = gamma,
                    beta0 = beta0,
                    beta1 = beta1,
                    mu = mu,
                    nu_eps = nu_eps,
                    nu_alpha = nu_alpha,
                    alpha0 = alpha0,
                    alpha1 = alpha1,
                    weights = weights))
}
