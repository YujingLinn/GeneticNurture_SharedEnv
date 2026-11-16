library(boot)
library(dplyr)
library(tidyr)
library(openxlsx)
library(parallel)
library(sandwich)
library(lmtest)

sourceFileStem <- "GN_Data/"
outFileStem <- "GN_Results/Reg/"
outFileStemRDS <- "GN_Results/Reg/RDS_output/"

if (!dir.exists(outFileStem))    dir.create(outFileStem,    recursive = TRUE)
if (!dir.exists(outFileStemRDS)) dir.create(outFileStemRDS, recursive = TRUE)

dat_scaled <- read.csv(paste(sourceFileStem, "GN_scaled.csv", sep=""))
source("0GN_VarPairs.R")

nboot <- 1000
ncpus <- 64

gn_boot_fn <- function(data, indices, outcome_var, child_pgs, mother_pgs, father_pgs, cov_str) {
  d <- data[indices, ]

  form_mat_pop <- as.formula(paste(outcome_var, "~", mother_pgs, "+", cov_str))
  form_pat_pop <- as.formula(paste(outcome_var, "~", father_pgs, "+", cov_str))
  form_covar_only  <- as.formula(paste(outcome_var, "~", cov_str))
  form_child_pop  <- as.formula(paste(outcome_var, "~", child_pgs, "+", cov_str))
  form_matGN <- as.formula(paste(outcome_var, "~", child_pgs, "+", mother_pgs, "+", cov_str))
  form_patGN <- as.formula(paste(outcome_var, "~", child_pgs, "+", father_pgs, "+", cov_str))
  form_overallGN  <- as.formula(paste(outcome_var, "~", child_pgs, "+", mother_pgs, "+", father_pgs, "+", cov_str))

  m_mat_pop <- lm(form_mat_pop, data = d)
  m_pat_pop <- lm(form_pat_pop, data = d)
  m_covar_only <- lm(form_covar_only, data = d)
  m_child_pop  <- lm(form_child_pop, data = d)
  m_matGN <- lm(form_matGN, data = d)
  m_patGN <- lm(form_patGN, data = d)
  m_overallGN  <- lm(form_overallGN, data = d)

  fix_mat_pop <- coef(m_mat_pop)
  fix_pat_pop <- coef(m_pat_pop)
  fix_child_pop  <- coef(m_child_pop)
  fix_matGN <- coef(m_matGN)
  fix_patGN <- coef(m_patGN)
  fix_overallGN  <- coef(m_overallGN)

  b_mom_pop        <- unname(fix_mat_pop[mother_pgs])
  b_dad_pop        <- unname(fix_pat_pop[father_pgs])
  b_child_pop      <- unname(fix_child_pop[child_pgs])
  b_child_matGN    <- unname(fix_matGN[child_pgs])
  b_mom_matGN      <- unname(fix_matGN[mother_pgs])
  b_child_patGN    <- unname(fix_patGN[child_pgs])
  b_dad_patGN      <- unname(fix_patGN[father_pgs])
  b_child_ovarllGN <- unname(fix_overallGN[child_pgs])
  b_mom_ovarllGN   <- unname(fix_overallGN[mother_pgs])
  b_dad_ovarllGN   <- unname(fix_overallGN[father_pgs])

  r2_cov_only <- summary(m_covar_only)$r.squared
  r2_mat_pop <- summary(m_mat_pop)$r.squared
  r2_pat_pop <- summary(m_pat_pop)$r.squared
  r2_child_pop  <- summary(m_child_pop)$r.squared
  r2_matGN <- summary(m_matGN)$r.squared
  r2_patGN <- summary(m_patGN)$r.squared
  r2_overallGN  <- summary(m_overallGN)$r.squared

  dr2_child <- r2_child_pop - r2_cov_only
  dr2_mat   <- r2_mat_pop  - r2_cov_only
  dr2_pat   <- r2_pat_pop  - r2_cov_only

  dr2_mat_marginal <- r2_matGN - r2_child_pop
  dr2_pat_marginal <- r2_patGN - r2_child_pop

  delta_r2_overall <- r2_overallGN - r2_child_pop
  delta_r2_uniq_mat <- r2_overallGN - r2_patGN
  delta_r2_uniq_pat <- r2_overallGN - r2_matGN
  delta_r2_shared_am <- delta_r2_overall - (delta_r2_uniq_mat + delta_r2_uniq_pat)

  denom_shrink <- b_child_pop
  denom_direct <- b_child_ovarllGN
  denom_share  <- b_mom_ovarllGN + b_dad_ovarllGN + b_child_ovarllGN

  shrinkage_child <- if (abs(denom_shrink) < 1e-8) NA_real_ else
    (b_child_pop - b_child_ovarllGN) / denom_shrink
  ratio_gn_direct <- if (abs(denom_direct) < 1e-8) NA_real_ else
    (b_mom_ovarllGN + b_dad_ovarllGN) / denom_direct
  prop_mat_gn     <- if (abs(denom_share)  < 1e-8) NA_real_ else
    b_mom_ovarllGN / denom_share
  prop_pat_gn     <- if (abs(denom_share)  < 1e-8) NA_real_ else
    b_dad_ovarllGN / denom_share

  c(
    Beta_Mom_Pop = b_mom_pop, Beta_Dad_Pop = b_dad_pop,
    Beta_Child_Pop = b_child_pop, Beta_Child_MatGN = b_child_matGN, Beta_Mom_MatGN = b_mom_matGN,
    Beta_Child_PatGN = b_child_patGN, Beta_Dad_PatGN = b_dad_patGN,
    Beta_Child_OverallGN = b_child_ovarllGN, Beta_Mom_OverallGN = b_mom_ovarllGN, Beta_Dad_OverallGN = b_dad_ovarllGN,
    R2_Covar_Only = r2_cov_only, R2_Mat_Pop = r2_mat_pop, R2_Pat_Pop = r2_pat_pop,
    R2_Child_Pop = r2_child_pop, R2_MatGN = r2_matGN, R2_PatGN = r2_patGN, R2_OverallGN = r2_overallGN,
    DeltaR2_Child = dr2_child, DeltaR2_Mat = dr2_mat, DeltaR2_Pat = dr2_pat,
    DeltaR2_Mat_Marginal = dr2_mat_marginal, DeltaR2_Pat_Marginal = dr2_pat_marginal,
    DeltaR2_Overall_GN = delta_r2_overall, DeltaR2_Unique_Mat = delta_r2_uniq_mat,
    DeltaR2_Unique_Pat = delta_r2_uniq_pat, DeltaR2_Shared_AM = delta_r2_shared_am,
    Shrinkage_Child = shrinkage_child, Ratio_GN_Direct = ratio_gn_direct,
    Prop_Mat_GN = prop_mat_gn, Prop_Pat_GN = prop_pat_gn
  )
}

run_gn_pipeline <- function(data, vp, child_pgs, mother_pgs, father_pgs, cov_str) {

  outcome_var <- paste0(vp[1], "1")
  n_used <- nrow(data)
  n_fam  <- length(unique(data$id_fam))

  s_mom <- sd(data[[mother_pgs]], na.rm = TRUE)
  s_dad <- sd(data[[father_pgs]], na.rm = TRUE)

  form_child_pop  <- as.formula(paste(outcome_var, "~", child_pgs, "+", cov_str))
  form_matGN <- as.formula(paste(outcome_var, "~", child_pgs, "+", mother_pgs, "+", cov_str))
  form_patGN <- as.formula(paste(outcome_var, "~", child_pgs, "+", father_pgs, "+", cov_str))
  form_overallGN  <- as.formula(paste(outcome_var, "~", child_pgs, "+", mother_pgs, "+", father_pgs, "+", cov_str))

  form_mat_pop <- as.formula(paste(outcome_var, "~", mother_pgs, "+", cov_str))
  form_pat_pop <- as.formula(paste(outcome_var, "~", father_pgs, "+", cov_str))
  form_covar_only <- as.formula(paste(outcome_var, "~", cov_str))

  m_mat_pop <- lm(form_mat_pop, data = data)
  m_pat_pop <- lm(form_pat_pop, data = data)
  m_covar_only <- lm(form_covar_only, data = data)
  m_child_pop  <- lm(form_child_pop, data = data)
  m_matGN <- lm(form_matGN, data = data)
  m_patGN <- lm(form_patGN, data = data)
  m_overallGN  <- lm(form_overallGN, data = data)

  cl <- data$id_fam
  crob <- function(m) lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = cl, type = "HC1"))

  sum_mat_pop <- crob(m_mat_pop)
  sum_pat_pop <- crob(m_pat_pop)
  sum_child_pop  <- crob(m_child_pop)
  sum_matGN <- crob(m_matGN)
  sum_patGN <- crob(m_patGN)
  sum_overallGN  <- crob(m_overallGN)

  se_mom_pop       <- sum_mat_pop[mother_pgs, "Std. Error"]
  se_dad_pop       <- sum_pat_pop[father_pgs, "Std. Error"]
  se_child_pop     <- sum_child_pop[child_pgs, "Std. Error"]
  se_child_matGN   <- sum_matGN[child_pgs, "Std. Error"]
  se_mom_matGN     <- sum_matGN[mother_pgs, "Std. Error"]
  se_child_patGN   <- sum_patGN[child_pgs, "Std. Error"]
  se_dad_patGN     <- sum_patGN[father_pgs, "Std. Error"]
  se_child_overall <- sum_overallGN[child_pgs, "Std. Error"]
  se_mom_overall   <- sum_overallGN[mother_pgs, "Std. Error"]
  se_dad_overall   <- sum_overallGN[father_pgs, "Std. Error"]

  p_mom_pop       <- sum_mat_pop[mother_pgs, "Pr(>|t|)"]
  p_dad_pop       <- sum_pat_pop[father_pgs, "Pr(>|t|)"]
  p_child_pop     <- sum_child_pop[child_pgs, "Pr(>|t|)"]
  p_child_matGN   <- sum_matGN[child_pgs, "Pr(>|t|)"]
  p_mom_matGN     <- sum_matGN[mother_pgs, "Pr(>|t|)"]
  p_child_patGN   <- sum_patGN[child_pgs, "Pr(>|t|)"]
  p_dad_patGN     <- sum_patGN[father_pgs, "Pr(>|t|)"]
  p_child_overall <- sum_overallGN[child_pgs, "Pr(>|t|)"]
  p_mom_overall   <- sum_overallGN[mother_pgs, "Pr(>|t|)"]
  p_dad_overall   <- sum_overallGN[father_pgs, "Pr(>|t|)"]

  vcov_overallGN   <- sandwich::vcovCL(m_overallGN, cluster = cl, type = "HC1")
  share_terms      <- c(mother_pgs, father_pgs, child_pgs)
  se_share_overall <- sqrt(sum(vcov_overallGN[share_terms, share_terms]))
  b_share_overall  <- sum_overallGN[mother_pgs, "Estimate"] + sum_overallGN[father_pgs, "Estimate"] + sum_overallGN[child_pgs, "Estimate"]
  p_share_overall  <- 2 * pt(-abs(b_share_overall / se_share_overall), df = m_overallGN$df.residual)

  vcl <- function(m) sandwich::vcovCL(m, cluster = cl, type = "HC1")
  waldf_child    <- lmtest::waldtest(m_covar_only, m_child_pop, vcov = vcl)$`Pr(>F)`[2]
  waldf_mat      <- lmtest::waldtest(m_covar_only, m_mat_pop,  vcov = vcl)$`Pr(>F)`[2]
  waldf_pat      <- lmtest::waldtest(m_covar_only, m_pat_pop,  vcov = vcl)$`Pr(>F)`[2]
  waldf_mat_marginal <- lmtest::waldtest(m_child_pop, m_matGN, vcov = vcl)$`Pr(>F)`[2]
  waldf_pat_marginal <- lmtest::waldtest(m_child_pop, m_patGN, vcov = vcl)$`Pr(>F)`[2]
  waldf_uniq_mat <- lmtest::waldtest(m_patGN, m_overallGN, vcov = vcl)$`Pr(>F)`[2]
  waldf_uniq_pat <- lmtest::waldtest(m_matGN, m_overallGN, vcov = vcl)$`Pr(>F)`[2]
  waldf_overall  <- lmtest::waldtest(m_child_pop, m_overallGN, vcov = vcl)$`Pr(>F)`[2]

  cf_ov  <- coef(m_overallGN)
  V_ov   <- vcl(m_overallGN)
  L_ov   <- setNames(rep(0, length(cf_ov)), names(cf_ov))
  L_ov[mother_pgs] <-  1
  L_ov[father_pgs] <- -1
  diff_mp_est <- sum(L_ov * cf_ov)
  diff_mp_se  <- sqrt(as.numeric(t(L_ov) %*% V_ov %*% L_ov))
  diff_mp_z   <- diff_mp_est / diff_mp_se
  diff_mp_p   <- 2 * pnorm(-abs(diff_mp_z))

  set.seed(2026)
  fam_ids  <- unique(data$id_fam)
  fam_rows <- split(seq_len(nrow(data)), data$id_fam)

  cluster_boot_fn <- function(fam_vec, idx, full_data, fam_rows,
                              outcome_var, child_pgs, mother_pgs, father_pgs, cov_str) {
    rows <- unlist(fam_rows[as.character(fam_vec[idx])], use.names = FALSE)
    gn_boot_fn(full_data, rows, outcome_var, child_pgs, mother_pgs, father_pgs, cov_str)
  }

  boot_out <- boot(data = fam_ids, statistic = cluster_boot_fn, R = nboot,
                   full_data = data, fam_rows = fam_rows,
                   outcome_var = outcome_var,
                   child_pgs = child_pgs, mother_pgs = mother_pgs, father_pgs = father_pgs, cov_str = cov_str,
                   parallel = "multicore", ncpus = ncpus)

  est <- boot_out$t0

  alpha_ratio_guard <- 0.05
  guard_ok <- c(
    Shrinkage_Child = (p_child_pop     <= alpha_ratio_guard),
    Ratio_GN_Direct = (p_child_overall <= alpha_ratio_guard),
    Prop_Mat_GN     = (p_share_overall <= alpha_ratio_guard),
    Prop_Pat_GN     = (p_share_overall <= alpha_ratio_guard)
  )

  ratio_names <- c("Shrinkage_Child", "Ratio_GN_Direct", "Prop_Mat_GN", "Prop_Pat_GN")
  for (nm in ratio_names) if (!guard_ok[[nm]]) est[[nm]] <- NA_real_
  ci_pb <- sapply(ratio_names, function(nm) {
    if (!guard_ok[[nm]]) return(c(NA_real_, NA_real_))
    boot.ci(boot_out, type = "perc", index = match(nm, names(est)))$percent[4:5]
  })
  colnames(ci_pb) <- ratio_names
  lo <- ci_pb[1, ]; hi <- ci_pb[2, ]
  z  <- qnorm(0.975)

  out_row <- data.frame(
    Base_Trait_Code = vp[1], PGS_Base = vp[2], Age_Var = vp[3],
    Trait_Category = vp[4], Domain = vp[5], Unique_Label = vp[6], Trait_Name = vp[7],
    Outcome_Var = outcome_var, N_used = n_used, N_fam = n_fam,

    R2_Covar_Only_Est = est[["R2_Covar_Only"]],

    Beta_Mom_Pop_Est = est[["Beta_Mom_Pop"]], Beta_Mom_Pop_SE = se_mom_pop,
    Beta_Mom_Pop_LowerCI = est[["Beta_Mom_Pop"]] - z*se_mom_pop, Beta_Mom_Pop_UpperCI = est[["Beta_Mom_Pop"]] + z*se_mom_pop,
    Beta_Mom_Pop_P = p_mom_pop,
    Beta_Mom_Pop_r_Est = est[["Beta_Mom_Pop"]] * s_mom, Beta_Mom_Pop_r_SE = se_mom_pop * s_mom,
    Beta_Mom_Pop_r_LowerCI = (est[["Beta_Mom_Pop"]] - z*se_mom_pop) * s_mom, Beta_Mom_Pop_r_UpperCI = (est[["Beta_Mom_Pop"]] + z*se_mom_pop) * s_mom,
    R2_Mat_Pop_Est = est[["R2_Mat_Pop"]],

    Beta_Dad_Pop_Est = est[["Beta_Dad_Pop"]], Beta_Dad_Pop_SE = se_dad_pop,
    Beta_Dad_Pop_LowerCI = est[["Beta_Dad_Pop"]] - z*se_dad_pop, Beta_Dad_Pop_UpperCI = est[["Beta_Dad_Pop"]] + z*se_dad_pop,
    Beta_Dad_Pop_P = p_dad_pop,
    Beta_Dad_Pop_r_Est = est[["Beta_Dad_Pop"]] * s_dad, Beta_Dad_Pop_r_SE = se_dad_pop * s_dad,
    Beta_Dad_Pop_r_LowerCI = (est[["Beta_Dad_Pop"]] - z*se_dad_pop) * s_dad, Beta_Dad_Pop_r_UpperCI = (est[["Beta_Dad_Pop"]] + z*se_dad_pop) * s_dad,
    R2_Pat_Pop_Est = est[["R2_Pat_Pop"]],

    Beta_Child_Pop_Est = est[["Beta_Child_Pop"]], Beta_Child_Pop_SE = se_child_pop,
    Beta_Child_Pop_LowerCI = est[["Beta_Child_Pop"]] - z*se_child_pop, Beta_Child_Pop_UpperCI = est[["Beta_Child_Pop"]] + z*se_child_pop,
    Beta_Child_Pop_P = p_child_pop,
    R2_Child_Pop_Est = est[["R2_Child_Pop"]],

    Beta_Child_MatGN_Est = est[["Beta_Child_MatGN"]], Beta_Child_MatGN_SE = se_child_matGN,
    Beta_Child_MatGN_LowerCI = est[["Beta_Child_MatGN"]] - z*se_child_matGN, Beta_Child_MatGN_UpperCI = est[["Beta_Child_MatGN"]] + z*se_child_matGN,
    Beta_Child_MatGN_P = p_child_matGN,
    Beta_Mom_MatGN_Est = est[["Beta_Mom_MatGN"]], Beta_Mom_MatGN_SE = se_mom_matGN,
    Beta_Mom_MatGN_LowerCI = est[["Beta_Mom_MatGN"]] - z*se_mom_matGN, Beta_Mom_MatGN_UpperCI = est[["Beta_Mom_MatGN"]] + z*se_mom_matGN,
    Beta_Mom_MatGN_P = p_mom_matGN,
    Beta_Mom_MatGN_r_Est = est[["Beta_Mom_MatGN"]] * s_mom, Beta_Mom_MatGN_r_SE = se_mom_matGN * s_mom,
    Beta_Mom_MatGN_r_LowerCI = (est[["Beta_Mom_MatGN"]] - z*se_mom_matGN) * s_mom, Beta_Mom_MatGN_r_UpperCI = (est[["Beta_Mom_MatGN"]] + z*se_mom_matGN) * s_mom,
    R2_MatGN_Est = est[["R2_MatGN"]],

    Beta_Child_PatGN_Est = est[["Beta_Child_PatGN"]], Beta_Child_PatGN_SE = se_child_patGN,
    Beta_Child_PatGN_LowerCI = est[["Beta_Child_PatGN"]] - z*se_child_patGN, Beta_Child_PatGN_UpperCI = est[["Beta_Child_PatGN"]] + z*se_child_patGN,
    Beta_Child_PatGN_P = p_child_patGN,
    Beta_Dad_PatGN_Est = est[["Beta_Dad_PatGN"]], Beta_Dad_PatGN_SE = se_dad_patGN,
    Beta_Dad_PatGN_LowerCI = est[["Beta_Dad_PatGN"]] - z*se_dad_patGN, Beta_Dad_PatGN_UpperCI = est[["Beta_Dad_PatGN"]] + z*se_dad_patGN,
    Beta_Dad_PatGN_P = p_dad_patGN,
    Beta_Dad_PatGN_r_Est = est[["Beta_Dad_PatGN"]] * s_dad, Beta_Dad_PatGN_r_SE = se_dad_patGN * s_dad,
    Beta_Dad_PatGN_r_LowerCI = (est[["Beta_Dad_PatGN"]] - z*se_dad_patGN) * s_dad, Beta_Dad_PatGN_r_UpperCI = (est[["Beta_Dad_PatGN"]] + z*se_dad_patGN) * s_dad,
    R2_PatGN_Est = est[["R2_PatGN"]],

    Beta_Child_OverallGN_Est = est[["Beta_Child_OverallGN"]], Beta_Child_OverallGN_SE = se_child_overall,
    Beta_Child_OverallGN_LowerCI = est[["Beta_Child_OverallGN"]] - z*se_child_overall, Beta_Child_OverallGN_UpperCI = est[["Beta_Child_OverallGN"]] + z*se_child_overall,
    Beta_Child_OverallGN_P = p_child_overall,
    Beta_Mom_OverallGN_Est = est[["Beta_Mom_OverallGN"]], Beta_Mom_OverallGN_SE = se_mom_overall,
    Beta_Mom_OverallGN_LowerCI = est[["Beta_Mom_OverallGN"]] - z*se_mom_overall, Beta_Mom_OverallGN_UpperCI = est[["Beta_Mom_OverallGN"]] + z*se_mom_overall,
    Beta_Mom_OverallGN_P = p_mom_overall,
    Beta_Mom_OverallGN_r_Est = est[["Beta_Mom_OverallGN"]] * s_mom, Beta_Mom_OverallGN_r_SE = se_mom_overall * s_mom,
    Beta_Mom_OverallGN_r_LowerCI = (est[["Beta_Mom_OverallGN"]] - z*se_mom_overall) * s_mom, Beta_Mom_OverallGN_r_UpperCI = (est[["Beta_Mom_OverallGN"]] + z*se_mom_overall) * s_mom,
    Beta_Dad_OverallGN_Est = est[["Beta_Dad_OverallGN"]], Beta_Dad_OverallGN_SE = se_dad_overall,
    Beta_Dad_OverallGN_LowerCI = est[["Beta_Dad_OverallGN"]] - z*se_dad_overall, Beta_Dad_OverallGN_UpperCI = est[["Beta_Dad_OverallGN"]] + z*se_dad_overall,
    Beta_Dad_OverallGN_P = p_dad_overall,
    Beta_Dad_OverallGN_r_Est = est[["Beta_Dad_OverallGN"]] * s_dad, Beta_Dad_OverallGN_r_SE = se_dad_overall * s_dad,
    Beta_Dad_OverallGN_r_LowerCI = (est[["Beta_Dad_OverallGN"]] - z*se_dad_overall) * s_dad, Beta_Dad_OverallGN_r_UpperCI = (est[["Beta_Dad_OverallGN"]] + z*se_dad_overall) * s_dad,
    R2_OverallGN_Est = est[["R2_OverallGN"]],

    DeltaR2_Child_Est = est[["DeltaR2_Child"]], WaldF_Child_P = waldf_child,
    DeltaR2_Mat_Est = est[["DeltaR2_Mat"]], WaldF_Mat_P = waldf_mat,
    DeltaR2_Pat_Est = est[["DeltaR2_Pat"]], WaldF_Pat_P = waldf_pat,
    DeltaR2_Mat_Marginal_Est = est[["DeltaR2_Mat_Marginal"]], WaldF_Mat_Marginal_P = waldf_mat_marginal,
    DeltaR2_Pat_Marginal_Est = est[["DeltaR2_Pat_Marginal"]], WaldF_Pat_Marginal_P = waldf_pat_marginal,
    DeltaR2_Unique_Mat_Est = est[["DeltaR2_Unique_Mat"]], WaldF_Unique_Mat_P = waldf_uniq_mat,
    DeltaR2_Unique_Pat_Est = est[["DeltaR2_Unique_Pat"]], WaldF_Unique_Pat_P = waldf_uniq_pat,
    DeltaR2_Overall_GN_Est = est[["DeltaR2_Overall_GN"]], WaldF_Overall_P = waldf_overall,
    DeltaR2_Shared_AM_Est = est[["DeltaR2_Shared_AM"]],

    Diff_Mom_Pat_Overall_Est = diff_mp_est, Diff_Mom_Pat_Overall_SE = diff_mp_se,
    Diff_Mom_Pat_Overall_P = diff_mp_p,

    Shrinkage_Child_Est = est[["Shrinkage_Child"]], Shrinkage_Child_LowerCI = lo[["Shrinkage_Child"]], Shrinkage_Child_UpperCI = hi[["Shrinkage_Child"]],
    Ratio_GN_Direct_Est = est[["Ratio_GN_Direct"]], Ratio_GN_Direct_LowerCI = lo[["Ratio_GN_Direct"]], Ratio_GN_Direct_UpperCI = hi[["Ratio_GN_Direct"]],
    Prop_Mat_GN_Est = est[["Prop_Mat_GN"]], Prop_Mat_GN_LowerCI = lo[["Prop_Mat_GN"]], Prop_Mat_GN_UpperCI = hi[["Prop_Mat_GN"]],
    Prop_Pat_GN_Est = est[["Prop_Pat_GN"]], Prop_Pat_GN_LowerCI = lo[["Prop_Pat_GN"]], Prop_Pat_GN_UpperCI = hi[["Prop_Pat_GN"]],
    stringsAsFactors = FALSE, check.names = FALSE
  )

  return(list(results_df = out_row, boot_object = boot_out))
}

run_duo_full_pipeline <- function(data, vp, child_pgs, parent_pgs, cov_str, parent_label) {

  outcome_var <- paste0(vp[1], "1")
  n_used <- nrow(data)
  n_fam  <- length(unique(data$id_fam))

  form_covar_only <- as.formula(paste(outcome_var, "~", cov_str))
  form_child_pop  <- as.formula(paste(outcome_var, "~", child_pgs, "+", cov_str))
  form_parent_pop <- as.formula(paste(outcome_var, "~", parent_pgs, "+", cov_str))
  form_duoGN      <- as.formula(paste(outcome_var, "~", child_pgs, "+", parent_pgs, "+", cov_str))

  m_covar_only <- lm(form_covar_only, data = data)
  m_child_pop  <- lm(form_child_pop,  data = data)
  m_parent_pop <- lm(form_parent_pop, data = data)
  m_duoGN      <- lm(form_duoGN,      data = data)

  cl    <- data$id_fam
  crob  <- function(m) lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = cl, type = "HC1"))
  vcl   <- function(m) sandwich::vcovCL(m, cluster = cl, type = "HC1")

  sum_child_pop  <- crob(m_child_pop)
  sum_parent_pop <- crob(m_parent_pop)
  sum_duoGN      <- crob(m_duoGN)
  z <- qnorm(0.975)

  s_parent <- sd(data[[parent_pgs]], na.rm = TRUE)

  beta_child_pop_est <- sum_child_pop[child_pgs, "Estimate"]
  beta_child_pop_se  <- sum_child_pop[child_pgs, "Std. Error"]
  beta_child_pop_p   <- sum_child_pop[child_pgs, "Pr(>|t|)"]

  beta_parent_pop_est <- sum_parent_pop[parent_pgs, "Estimate"]
  beta_parent_pop_se  <- sum_parent_pop[parent_pgs, "Std. Error"]
  beta_parent_pop_p   <- sum_parent_pop[parent_pgs, "Pr(>|t|)"]

  beta_child_duo_est <- sum_duoGN[child_pgs, "Estimate"]
  beta_child_duo_se  <- sum_duoGN[child_pgs, "Std. Error"]
  beta_child_duo_p   <- sum_duoGN[child_pgs, "Pr(>|t|)"]

  beta_parent_est <- sum_duoGN[parent_pgs, "Estimate"]
  beta_parent_se  <- sum_duoGN[parent_pgs, "Std. Error"]
  beta_parent_p   <- sum_duoGN[parent_pgs, "Pr(>|t|)"]

  waldf_child <- lmtest::waldtest(m_covar_only, m_child_pop, vcov = vcl)$`Pr(>F)`[2]
  waldf_parent_pop <- lmtest::waldtest(m_covar_only, m_parent_pop, vcov = vcl)$`Pr(>F)`[2]
  waldf_marginal <- lmtest::waldtest(m_child_pop, m_duoGN, vcov = vcl)$`Pr(>F)`[2]

  r2_covar_only <- summary(m_covar_only)$r.squared
  r2_child_pop  <- summary(m_child_pop)$r.squared
  r2_parent_pop <- summary(m_parent_pop)$r.squared
  r2_duoGN      <- summary(m_duoGN)$r.squared

  gn_tag      <- if (parent_label == "Mom") "MatGN" else "PatGN"
  marg_tag    <- if (parent_label == "Mom") "Mat_Marginal" else "Pat_Marginal"
  pop_r2_tag  <- if (parent_label == "Mom") "Mat_Pop" else "Pat_Pop"
  pop_delta_tag <- if (parent_label == "Mom") "Mat" else "Pat"
  parent_pop_name  <- paste0("Beta_", parent_label, "_Pop")
  parent_duo_name  <- paste0("Beta_", parent_label, "_", gn_tag)

  out_row <- data.frame(
    Outcome_Var = outcome_var,
    N_used = n_used,
    N_fam  = n_fam,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  out_row[["R2_Covar_Only_Est"]]      <- r2_covar_only
  out_row[["Beta_Child_Pop_Est"]]     <- beta_child_pop_est
  out_row[["Beta_Child_Pop_SE"]]      <- beta_child_pop_se
  out_row[["Beta_Child_Pop_LowerCI"]] <- beta_child_pop_est - z * beta_child_pop_se
  out_row[["Beta_Child_Pop_UpperCI"]] <- beta_child_pop_est + z * beta_child_pop_se
  out_row[["Beta_Child_Pop_P"]]       <- beta_child_pop_p
  out_row[["R2_Child_Pop_Est"]]       <- r2_child_pop
  out_row[["DeltaR2_Child_Est"]]      <- r2_child_pop - r2_covar_only
  out_row[["WaldF_Child_P"]]          <- waldf_child

  out_row[[paste0(parent_pop_name, "_Est")]]     <- beta_parent_pop_est
  out_row[[paste0(parent_pop_name, "_SE")]]      <- beta_parent_pop_se
  out_row[[paste0(parent_pop_name, "_LowerCI")]] <- beta_parent_pop_est - z * beta_parent_pop_se
  out_row[[paste0(parent_pop_name, "_UpperCI")]] <- beta_parent_pop_est + z * beta_parent_pop_se
  out_row[[paste0(parent_pop_name, "_P")]]       <- beta_parent_pop_p
  out_row[[paste0(parent_pop_name, "_r_Est")]]     <- beta_parent_pop_est * s_parent
  out_row[[paste0(parent_pop_name, "_r_SE")]]      <- beta_parent_pop_se * s_parent
  out_row[[paste0(parent_pop_name, "_r_LowerCI")]] <- (beta_parent_pop_est - z * beta_parent_pop_se) * s_parent
  out_row[[paste0(parent_pop_name, "_r_UpperCI")]] <- (beta_parent_pop_est + z * beta_parent_pop_se) * s_parent
  out_row[[paste0("R2_", pop_r2_tag, "_Est")]]        <- r2_parent_pop
  out_row[[paste0("DeltaR2_", pop_delta_tag, "_Est")]] <- r2_parent_pop - r2_covar_only
  out_row[[paste0("WaldF_", pop_delta_tag, "_P")]]     <- waldf_parent_pop

  out_row[[paste0("Beta_Child_", gn_tag, "_Est")]]     <- beta_child_duo_est
  out_row[[paste0("Beta_Child_", gn_tag, "_SE")]]      <- beta_child_duo_se
  out_row[[paste0("Beta_Child_", gn_tag, "_LowerCI")]] <- beta_child_duo_est - z * beta_child_duo_se
  out_row[[paste0("Beta_Child_", gn_tag, "_UpperCI")]] <- beta_child_duo_est + z * beta_child_duo_se
  out_row[[paste0("Beta_Child_", gn_tag, "_P")]]       <- beta_child_duo_p

  out_row[[paste0(parent_duo_name, "_Est")]]     <- beta_parent_est
  out_row[[paste0(parent_duo_name, "_SE")]]      <- beta_parent_se
  out_row[[paste0(parent_duo_name, "_LowerCI")]] <- beta_parent_est - z * beta_parent_se
  out_row[[paste0(parent_duo_name, "_UpperCI")]] <- beta_parent_est + z * beta_parent_se
  out_row[[paste0(parent_duo_name, "_P")]]       <- beta_parent_p
  out_row[[paste0(parent_duo_name, "_r_Est")]]     <- beta_parent_est * s_parent
  out_row[[paste0(parent_duo_name, "_r_SE")]]      <- beta_parent_se * s_parent
  out_row[[paste0(parent_duo_name, "_r_LowerCI")]] <- (beta_parent_est - z * beta_parent_se) * s_parent
  out_row[[paste0(parent_duo_name, "_r_UpperCI")]] <- (beta_parent_est + z * beta_parent_se) * s_parent

  out_row[[paste0("R2_", gn_tag, "_Est")]]        <- r2_duoGN
  out_row[[paste0("DeltaR2_", marg_tag, "_Est")]] <- r2_duoGN - r2_child_pop
  out_row[[paste0("WaldF_", marg_tag, "_P")]]     <- waldf_marginal

  return(out_row)
}

run_subgroup_analysis <- function(data, var_pairs_list, subgroup_label = "All_Trios", include_sex = TRUE, pgs_suffix = "") {

  all_results <- list()
  all_boot_objs <- list()

  for (vp in var_pairs_list) {
    outcome_var <- paste0(vp[1], "1")
    pgs_base    <- vp[2]
    age_var     <- vp[3]
    unique_lbl  <- vp[6]

    cat(sprintf("Running [%s] for %s...\n", subgroup_label, unique_lbl))

    child_pgs  <- paste0(pgs_base, pgs_suffix)
    mother_pgs <- paste0(pgs_base, "_fp", pgs_suffix)
    father_pgs <- paste0(pgs_base, "_mp", pgs_suffix)

    covariates <- c(age_var, "chiptype")
    if (include_sex == TRUE) { covariates <- c(covariates, "sex1") }
    cov_str <- paste(ifelse(covariates == "chiptype", "as.factor(chiptype)", covariates), collapse = " + ")

    vars_needed <- c(outcome_var, child_pgs, mother_pgs, father_pgs, covariates, "id_fam")
    missing_cols <- setdiff(vars_needed, names(data))

    if(length(missing_cols) > 0) {
      cat(sprintf("-> SKIPPING: Missing columns %s\n\n", paste(missing_cols, collapse = ", ")))
      next
    }

    df_clean <- data[complete.cases(data[, vars_needed]), ]

    if(nrow(df_clean) < 100) {
      cat(sprintf("-> SKIPPING due to low sample size (N=%d)\n\n", nrow(df_clean)))
      next
    }

    pipeline_out <- run_gn_pipeline(data = df_clean, vp = vp,
                                    child_pgs = child_pgs, mother_pgs = mother_pgs, father_pgs = father_pgs,
                                    cov_str = cov_str)

    res_df <- pipeline_out$results_df %>%
      mutate(Subgroup = subgroup_label) %>%
      relocate(Subgroup, .after = Outcome_Var)
    all_results[[unique_lbl]] <- res_df

    all_boot_objs[[unique_lbl]] <- pipeline_out$boot_object

    cat("-> Success!\n\n")
  }

  rds_filename <- paste0(outFileStemRDS, "boot_objects_", subgroup_label, ".rds")
  if(!dir.exists(outFileStemRDS)) dir.create(outFileStemRDS, recursive = TRUE)
  saveRDS(all_boot_objs, file = rds_filename)
  cat(sprintf("-> Saved RDS file for %s to %s\n", subgroup_label, rds_filename))

  return(bind_rows(all_results))
}

run_subgroup_analysis_duo_full <- function(data, var_pairs_list, parent_side = c("Mom", "Dad"), include_sex = TRUE) {

  parent_side <- match.arg(parent_side)
  subgroup_label  <- if (parent_side == "Mom") "Maternal_Duo_Full" else "Paternal_Duo_Full"
  parent_suffix   <- if (parent_side == "Mom") "_fp" else "_mp"

  all_results <- list()

  for (vp in var_pairs_list) {
    outcome_var <- paste0(vp[1], "1")
    pgs_base    <- vp[2]
    age_var     <- vp[3]
    unique_lbl  <- vp[6]

    cat(sprintf("Running [%s] for %s...\n", subgroup_label, unique_lbl))

    child_pgs  <- pgs_base
    parent_pgs <- paste0(pgs_base, parent_suffix)

    covariates <- c(age_var, "chiptype")
    if (include_sex == TRUE) { covariates <- c(covariates, "sex1") }
    cov_str <- paste(ifelse(covariates == "chiptype", "as.factor(chiptype)", covariates), collapse = " + ")

    vars_needed <- c(outcome_var, child_pgs, parent_pgs, covariates, "id_fam")
    missing_cols <- setdiff(vars_needed, names(data))

    if (length(missing_cols) > 0) {
      cat(sprintf("-> SKIPPING: Missing columns %s\n\n", paste(missing_cols, collapse = ", ")))
      next
    }

    df_clean <- data[complete.cases(data[, vars_needed]), ]

    if (nrow(df_clean) < 100) {
      cat(sprintf("-> SKIPPING due to low sample size (N=%d)\n\n", nrow(df_clean)))
      next
    }

    res_df <- run_duo_full_pipeline(data = df_clean, vp = vp,
                                    child_pgs = child_pgs, parent_pgs = parent_pgs,
                                    cov_str = cov_str, parent_label = parent_side) %>%
      mutate(Domain = vp[5], Unique_Label = unique_lbl, Subgroup = subgroup_label) %>%
      relocate(Domain, Unique_Label, Subgroup, .after = Outcome_Var)

    all_results[[unique_lbl]] <- res_df
    cat("-> Success!\n\n")
  }

  return(bind_rows(all_results))
}

run_subgroup_analysis_sibIGE <- function(data_dz, var_pairs_list, include_sex = TRUE) {

  subgroup_label <- "DZ_SibIGE_Quads"
  all_results <- list()
  all_boot_objs <- list()

  for (vp in var_pairs_list) {
    outcome_var <- paste0(vp[1], "1")
    pgs_base    <- vp[2]
    age_var     <- vp[3]
    unique_lbl  <- vp[6]

    cat(sprintf("Running [%s] for %s...\n", subgroup_label, unique_lbl))

    child_pgs  <- pgs_base
    mother_pgs <- paste0(pgs_base, "_fp")
    father_pgs <- paste0(pgs_base, "_mp")
    cotwin_pgs <- paste0(pgs_base, "_cotwin")

    covariates <- c(age_var, "chiptype")
    if (include_sex == TRUE) { covariates <- c(covariates, "sex1") }
    cov_str <- paste(ifelse(covariates == "chiptype", "as.factor(chiptype)", covariates), collapse = " + ")

    lk <- data_dz[, c("id_fam", "id_twin", child_pgs)]
    names(lk) <- c("id_fam", "cotwin_id", cotwin_pgs)
    d_sib <- merge(data_dz, lk, by = "id_fam")
    d_sib <- d_sib[d_sib$id_twin != d_sib$cotwin_id, ]

    vars_needed <- c(outcome_var, child_pgs, mother_pgs, father_pgs, cotwin_pgs, covariates, "id_fam")
    missing_cols <- setdiff(vars_needed, names(d_sib))

    if (length(missing_cols) > 0) {
      cat(sprintf("-> SKIPPING: Missing columns %s\n\n", paste(missing_cols, collapse = ", ")))
      next
    }

    df_clean <- d_sib[complete.cases(d_sib[, vars_needed]), ]

    if (nrow(df_clean) < 100) {
      cat(sprintf("-> SKIPPING due to low sample size (N=%d)\n\n", nrow(df_clean)))
      next
    }

    pipeline_out <- run_gn_pipeline(data = df_clean, vp = vp,
                                    child_pgs = child_pgs, mother_pgs = mother_pgs,
                                    father_pgs = father_pgs, cov_str = cov_str)

    form_sibGN <- as.formula(paste(outcome_var, "~", child_pgs, "+", mother_pgs, "+", father_pgs, "+", cotwin_pgs, "+", cov_str))
    m_sibGN <- lm(form_sibGN, data = df_clean)
    cl <- df_clean$id_fam
    sum_sibGN <- lmtest::coeftest(m_sibGN, vcov. = sandwich::vcovCL(m_sibGN, cluster = cl, type = "HC1"))
    z <- qnorm(0.975)

    s_mom <- sd(df_clean[[mother_pgs]], na.rm = TRUE)
    s_dad <- sd(df_clean[[father_pgs]], na.rm = TRUE)

    beta_child_sib  <- sum_sibGN[child_pgs,  "Estimate"]; se_child_sib  <- sum_sibGN[child_pgs,  "Std. Error"]; p_child_sib  <- sum_sibGN[child_pgs,  "Pr(>|t|)"]
    beta_mom_sib    <- sum_sibGN[mother_pgs, "Estimate"]; se_mom_sib    <- sum_sibGN[mother_pgs, "Std. Error"]; p_mom_sib    <- sum_sibGN[mother_pgs, "Pr(>|t|)"]
    beta_dad_sib    <- sum_sibGN[father_pgs, "Estimate"]; se_dad_sib    <- sum_sibGN[father_pgs, "Std. Error"]; p_dad_sib    <- sum_sibGN[father_pgs, "Pr(>|t|)"]
    beta_cotwin_sib <- sum_sibGN[cotwin_pgs, "Estimate"]; se_cotwin_sib <- sum_sibGN[cotwin_pgs, "Std. Error"]; p_cotwin_sib <- sum_sibGN[cotwin_pgs, "Pr(>|t|)"]
    r2_sibGN <- summary(m_sibGN)$r.squared

    sib_extra <- data.frame(
      Beta_Child_SibGN_Est = beta_child_sib, Beta_Child_SibGN_SE = se_child_sib,
      Beta_Child_SibGN_LowerCI = beta_child_sib - z*se_child_sib, Beta_Child_SibGN_UpperCI = beta_child_sib + z*se_child_sib,
      Beta_Child_SibGN_P = p_child_sib,

      Beta_Mom_SibGN_Est = beta_mom_sib, Beta_Mom_SibGN_SE = se_mom_sib,
      Beta_Mom_SibGN_LowerCI = beta_mom_sib - z*se_mom_sib, Beta_Mom_SibGN_UpperCI = beta_mom_sib + z*se_mom_sib,
      Beta_Mom_SibGN_P = p_mom_sib,
      Beta_Mom_SibGN_r_Est = beta_mom_sib * s_mom, Beta_Mom_SibGN_r_SE = se_mom_sib * s_mom,
      Beta_Mom_SibGN_r_LowerCI = (beta_mom_sib - z*se_mom_sib) * s_mom, Beta_Mom_SibGN_r_UpperCI = (beta_mom_sib + z*se_mom_sib) * s_mom,

      Beta_Dad_SibGN_Est = beta_dad_sib, Beta_Dad_SibGN_SE = se_dad_sib,
      Beta_Dad_SibGN_LowerCI = beta_dad_sib - z*se_dad_sib, Beta_Dad_SibGN_UpperCI = beta_dad_sib + z*se_dad_sib,
      Beta_Dad_SibGN_P = p_dad_sib,
      Beta_Dad_SibGN_r_Est = beta_dad_sib * s_dad, Beta_Dad_SibGN_r_SE = se_dad_sib * s_dad,
      Beta_Dad_SibGN_r_LowerCI = (beta_dad_sib - z*se_dad_sib) * s_dad, Beta_Dad_SibGN_r_UpperCI = (beta_dad_sib + z*se_dad_sib) * s_dad,

      Beta_Cotwin_SibGN_Est = beta_cotwin_sib, Beta_Cotwin_SibGN_SE = se_cotwin_sib,
      Beta_Cotwin_SibGN_LowerCI = beta_cotwin_sib - z*se_cotwin_sib, Beta_Cotwin_SibGN_UpperCI = beta_cotwin_sib + z*se_cotwin_sib,
      Beta_Cotwin_SibGN_P = p_cotwin_sib,

      R2_SibGN_Est = r2_sibGN,
      stringsAsFactors = FALSE, check.names = FALSE
    )

    res_df <- pipeline_out$results_df %>%
      mutate(Subgroup = subgroup_label) %>%
      relocate(Subgroup, .after = Outcome_Var) %>%
      bind_cols(sib_extra)

    all_results[[unique_lbl]]   <- res_df
    all_boot_objs[[unique_lbl]] <- pipeline_out$boot_object
    cat("-> Success!\n\n")
  }

  rds_filename <- paste0(outFileStemRDS, "boot_objects_", subgroup_label, ".rds")
  if(!dir.exists(outFileStemRDS)) dir.create(outFileStemRDS, recursive = TRUE)
  saveRDS(all_boot_objs, file = rds_filename)
  cat(sprintf("-> Saved RDS file for %s to %s\n", subgroup_label, rds_filename))

  return(bind_rows(all_results))
}

dat_female <- dat_scaled %>% filter(sex1 == 0)
dat_male   <- dat_scaled %>% filter(sex1 == 1)
dat_mz     <- dat_scaled %>% filter(zygos == 1)
dat_dz     <- dat_scaled %>% filter(zygos == 2)

dat_unrelated <- read.csv(paste(sourceFileStem, "GN_scaled_selectunpaired.csv", sep = ""))
res_whole  <- run_subgroup_analysis(dat_scaled, VarPairs, "All_Trios", include_sex = TRUE)
res_female <- run_subgroup_analysis(dat_female, VarPairs, "Female_Trios", include_sex = FALSE)
res_male   <- run_subgroup_analysis(dat_male,   VarPairs, "Male_Trios",   include_sex = FALSE)
res_mz     <- run_subgroup_analysis(dat_mz,     VarPairs, "MZ_Trios",     include_sex = TRUE)
res_dz     <- run_subgroup_analysis(dat_dz,     VarPairs, "DZ_Trios",     include_sex = TRUE)
res_unrelated <- run_subgroup_analysis(dat_unrelated, VarPairs, "Unrelated_Trios", include_sex = TRUE)
res_noPC   <- run_subgroup_analysis(dat_scaled, VarPairs, "noPC_Trios", include_sex = TRUE, pgs_suffix = "_noPC")

res_duo_full_mom <- run_subgroup_analysis_duo_full(dat_scaled, VarPairs, parent_side = "Mom", include_sex = TRUE)
res_duo_full_dad <- run_subgroup_analysis_duo_full(dat_scaled, VarPairs, parent_side = "Dad", include_sex = TRUE)

res_sibIGE <- run_subgroup_analysis_sibIGE(dat_dz, VarPairs, include_sex = TRUE)

subgroup_wide <- list(
  All_Trios          = res_whole,
  Female_Trios       = res_female,
  Male_Trios         = res_male,
  MZ_Trios           = res_mz,
  DZ_Trios           = res_dz,
  Unrelated_Trios    = res_unrelated,
  noPC_Trios         = res_noPC,
  Maternal_Duo_Full  = res_duo_full_mom,
  Paternal_Duo_Full  = res_duo_full_dad,
  DZ_SibIGE_Quads    = res_sibIGE
)

add_fdr <- function(w) {
  beta_p <- grep("^Beta_.*_P$", names(w), value = TRUE)
  stopifnot("Domain" %in% names(w))
  w <- w %>%
    group_by(Domain) %>%
    mutate(across(all_of(beta_p), ~ p.adjust(.x, method = "fdr"), .names = "{.col}_FDR")) %>%
    ungroup()
  for (p in beta_p) {
    w <- relocate(w, all_of(paste0(p, "_FDR")), .after = all_of(p))
  }
  w
}

wide_list <- lapply(subgroup_wide, add_fdr)
final_results_wide <- bind_rows(wide_list)

id_cols <- c("Domain", "Unique_Label", "Subgroup", "N_used")

table_params <- final_results_wide %>%
  select(all_of(id_cols),
         N_fam,
         starts_with("Beta_"),
         R2_Covar_Only_Est, R2_Mat_Pop_Est, R2_Pat_Pop_Est, R2_Child_Pop_Est,
         R2_MatGN_Est, R2_PatGN_Est, R2_OverallGN_Est, R2_SibGN_Est)

table_derived <- final_results_wide %>%
  select(all_of(id_cols),
         DeltaR2_Child_Est, WaldF_Child_P, DeltaR2_Mat_Est, WaldF_Mat_P, DeltaR2_Pat_Est, WaldF_Pat_P,
         DeltaR2_Mat_Marginal_Est, WaldF_Mat_Marginal_P, DeltaR2_Pat_Marginal_Est, WaldF_Pat_Marginal_P,
         DeltaR2_Unique_Mat_Est, WaldF_Unique_Mat_P, DeltaR2_Unique_Pat_Est, WaldF_Unique_Pat_P,
         DeltaR2_Overall_GN_Est, WaldF_Overall_P, DeltaR2_Shared_AM_Est,
         Diff_Mom_Pat_Overall_Est, Diff_Mom_Pat_Overall_SE, Diff_Mom_Pat_Overall_P,
         starts_with("Shrinkage_Child"), starts_with("Ratio_GN_Direct"),
         starts_with("Prop_Mat_GN"), starts_with("Prop_Pat_GN"))

table_all <- table_params %>% left_join(table_derived, by = id_cols)
write.csv(table_all, paste0(outFileStem, "GN_Reg.csv"), row.names = FALSE)

parent_beta_bases <- c("Beta_Mom_Pop", "Beta_Dad_Pop",
                       "Beta_Mom_MatGN", "Beta_Dad_PatGN",
                       "Beta_Mom_OverallGN", "Beta_Dad_OverallGN",
                       "Beta_Mom_SibGN", "Beta_Dad_SibGN")
excel_drop_cols <- as.vector(outer(parent_beta_bases,
                                   c("_Est", "_SE", "_LowerCI", "_UpperCI"),
                                   paste0))
table_params_excel <- table_params %>% select(-any_of(excel_drop_cols))

table_duo_full_excel <- table_all %>%
  filter(Subgroup %in% c("Maternal_Duo_Full", "Paternal_Duo_Full")) %>%
  select(where(~ !all(is.na(.))))

table_sibIGE_excel <- table_all %>%
  filter(Subgroup == "DZ_SibIGE_Quads") %>%
  select(where(~ !all(is.na(.))))

wb <- createWorkbook()
addWorksheet(wb, "Params");  writeData(wb, "Params",  table_params_excel);  freezePane(wb, "Params",  firstRow = TRUE)
addWorksheet(wb, "Derived"); writeData(wb, "Derived", table_derived);       freezePane(wb, "Derived", firstRow = TRUE)
addWorksheet(wb, "Duo_Full_Sensitivity"); writeData(wb, "Duo_Full_Sensitivity", table_duo_full_excel); freezePane(wb, "Duo_Full_Sensitivity", firstRow = TRUE)
addWorksheet(wb, "SibIGE_Quads"); writeData(wb, "SibIGE_Quads", table_sibIGE_excel); freezePane(wb, "SibIGE_Quads", firstRow = TRUE)
saveWorkbook(wb, paste0(outFileStem, "GN_Reg.xlsx"), overwrite = TRUE)

cat("Pipeline complete! Check your output folder.\n")
