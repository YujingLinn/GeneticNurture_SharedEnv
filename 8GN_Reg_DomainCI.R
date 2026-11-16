library(boot)
library(dplyr)
library(parallel)

sourceFileStem <- "GN_Data/"
outFileStem    <- "GN_Results/Reg/"
outFileStemRDS <- "GN_Results/Reg/RDS_output/"

dat_scaled <- read.csv(paste(sourceFileStem, "GN_scaled.csv", sep = ""))
source("0GN_VarPairs.R")

nboot <- 1000
ncpus <- 64

trait_list <- list()
for (vp in VarPairs) {
  outcome_var <- paste0(vp[1], "1")
  child_pgs   <- vp[2]
  mother_pgs  <- paste0(vp[2], "_fp")
  father_pgs  <- paste0(vp[2], "_mp")
  covariates  <- c(vp[3], "chiptype", "sex1")
  cov_str     <- paste(ifelse(covariates == "chiptype", "as.factor(chiptype)", covariates), collapse = " + ")

  vars_needed  <- c(outcome_var, child_pgs, mother_pgs, father_pgs, covariates, "id_fam")
  missing_cols <- setdiff(vars_needed, names(dat_scaled))
  if (length(missing_cols) > 0) { cat(sprintf("-> SKIPPING %s: missing %s\n", vp[6], paste(missing_cols, collapse = ", "))); next }

  df_clean <- dat_scaled[complete.cases(dat_scaled[, vars_needed]), ]
  if (nrow(df_clean) < 100) { cat(sprintf("-> SKIPPING %s (N=%d)\n", vp[6], nrow(df_clean))); next }

  trait_list[[vp[6]]] <- list(
    label      = vp[6],
    domain     = vp[7],
    data       = df_clean,
    rows_by_fam = split(seq_len(nrow(df_clean)), as.character(df_clean$id_fam)),
    child_pgs  = child_pgs, mother_pgs = mother_pgs, father_pgs = father_pgs,
    form_pop   = as.formula(paste(outcome_var, "~", child_pgs, "+", cov_str)),
    form_trio  = as.formula(paste(outcome_var, "~", child_pgs, "+", mother_pgs, "+", father_pgs, "+", cov_str)),
    s_mom      = sd(df_clean[[mother_pgs]], na.rm = TRUE),
    s_dad      = sd(df_clean[[father_pgs]], na.rm = TRUE)
  )
}
cat(sprintf("%d traits prepared\n", length(trait_list)))

fit_trait <- function(tr, d) {
  b_pop  <- coef(lm(tr$form_pop,  data = d))
  b_trio <- coef(lm(tr$form_trio, data = d))
  c(Pop    = unname(b_pop[tr$child_pgs]),
    Direct = unname(b_trio[tr$child_pgs]),
    Mom    = unname(b_trio[tr$mother_pgs]) * tr$s_mom,
    Dad    = unname(b_trio[tr$father_pgs]) * tr$s_dad)
}

qty <- c("Mom", "Dad", "Direct", "Pop")
domains <- unique(sapply(trait_list, `[[`, "domain"))
boot_objs  <- list()
results    <- list()

for (dm in domains) {
  trs <- trait_list[sapply(trait_list, `[[`, "domain") == dm]
  fam_union <- sort(unique(unlist(lapply(trs, function(tr) names(tr$rows_by_fam)))))
  cat(sprintf("Running domain [%s]: %d traits, %d families in union...\n", dm, length(trs), length(fam_union)))

  domain_stat <- function(fams, idx) {
    fam_b <- fams[idx]
    per_trait <- sapply(trs, function(tr) {
      rows <- unlist(tr$rows_by_fam[fam_b], use.names = FALSE)
      fit_trait(tr, tr$data[rows, ])
    })
    trait_vec <- as.vector(per_trait)
    names(trait_vec) <- paste(rep(colnames(per_trait), each = 4), rownames(per_trait), sep = "::")
    dom_mean <- rowMeans(per_trait)[qty]
    names(dom_mean) <- paste0("DomainMean_", qty)
    c(dom_mean, trait_vec)
  }

  set.seed(2026)
  b <- boot(data = fam_union, statistic = domain_stat, R = nboot,
            parallel = "multicore", ncpus = ncpus)
  boot_objs[[dm]] <- b

  ci <- sapply(seq_along(qty), function(j) boot.ci(b, type = "perc", index = j)$percent[4:5])
  results[[dm]] <- data.frame(
    Domain   = dm,
    N_traits = length(trs),
    N_fam_union = length(fam_union),
    Quantity = qty,
    Est      = unname(b$t0[seq_along(qty)]),
    LowerCI  = ci[1, ],
    UpperCI  = ci[2, ]
  )
}

domain_ci <- bind_rows(results) %>%
  mutate(Formatted = sprintf("%.2f (%.2f, %.2f)", Est, LowerCI, UpperCI))

table2_ci <- domain_ci %>%
  select(Domain, Quantity, Formatted) %>%
  tidyr::pivot_wider(names_from = Quantity, values_from = Formatted) %>%
  select(Domain, Mom, Dad, Direct, Pop)

write.csv(domain_ci, paste0(outFileStem, "GN_Table2_DomainCI.csv"), row.names = FALSE)
write.csv(table2_ci, paste0(outFileStem, "GN_Table2_DomainCI_formatted.csv"), row.names = FALSE)
saveRDS(list(boot_objs = boot_objs, nboot = nboot, session = sessionInfo()),
        paste0(outFileStemRDS, "GN_Table2_DomainCI_boot.rds"))
