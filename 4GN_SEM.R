library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)

sourceFileStem <- "GN_Data/"
outFileStem    <- "GN_Results/SEM/"
regFileStem    <- "GN_Results/Reg/"

if (!dir.exists(outFileStem)) dir.create(outFileStem, recursive = TRUE)

df_sem <- read.csv(paste0(sourceFileStem, "GN_scaled.csv"))
source("0GN_VarPairs.R")
source("0GN_Theme.R")

extract_param <- function(fit, param_labels, trait_label, model_name) {
  pe <- parameterEstimates(fit, ci = TRUE)
  pe <- pe[pe$label %in% param_labels, , drop = FALSE]

  spread <- tapply(pe$est, pe$label, function(x) diff(range(x, na.rm = TRUE)))
  if (any(spread > 1e-8, na.rm = TRUE))
    stop(sprintf("[%s | %s] Estimates for the same label differ across groups (max spread = %.3g); extract_param() must extract by group.",
                 model_name, trait_label, max(spread)))

  pe <- pe[!duplicated(pe$label), , drop = FALSE]

  results <- lapply(param_labels, function(lbl) {
    row <- pe[pe$label == lbl, ]
    if (nrow(row) == 0) return(NULL)
    data.frame(
      trait    = trait_label,
      model    = model_name,
      param    = lbl,
      est      = row$est,
      se       = row$se,
      pvalue   = row$pvalue,
      ci_lower = row$ci.lower,
      ci_upper = row$ci.upper,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(results)
}

extract_fit <- function(fit, trait_label, model_name) {
  fi <- fitMeasures(fit, c("chisq.scaled", "df", "pvalue.scaled",
                           "cfi.scaled", "rmsea.scaled",
                           "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled",
                           "cfi",
                           "srmr", "aic", "bic"))
  data.frame(
    trait    = trait_label,
    model    = model_name,
    converged = lavInspect(fit, "converged"),
    n_used    = lavInspect(fit, "ntotal"),
    chisq    = fi["chisq.scaled"],
    df       = fi["df"],
    p_chisq  = fi["pvalue.scaled"],
    cfi      = fi["cfi.scaled"],
    cfi_raw  = fi["cfi"],
    rmsea    = fi["rmsea.scaled"],
    rmsea_ll = fi["rmsea.ci.lower.scaled"],
    rmsea_ul = fi["rmsea.ci.upper.scaled"],
    srmr     = fi["srmr"],
    aic      = fi["aic"],
    bic      = fi["bic"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

build_df_mg <- function(df_sem, y_var, child_pgs, mom_pgs, dad_pgs, covs) {
  df_sem %>%
    filter(
      !is.na(.data[[child_pgs]]),
      !is.na(.data[[y_var]]),
      if_all(all_of(covs), ~ !is.na(.))
    ) %>%
    mutate(
      group = case_when(
        !is.na(.data[[mom_pgs]]) & !is.na(.data[[dad_pgs]]) ~ "trios",
        !is.na(.data[[mom_pgs]]) &  is.na(.data[[dad_pgs]]) ~ "dad_missing",
        is.na(.data[[mom_pgs]]) & !is.na(.data[[dad_pgs]]) ~ "mom_missing",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group))
}

all_results <- list()
fit_results <- list()

fit_one <- function(model_str, df_mg) {
  multi_group <- length(unique(df_mg$group)) >= 2
  has_cluster <- all(vapply(split(df_mg$id_fam, df_mg$group),
                            function(x) max(table(x)) > 1, logical(1)))
  args <- list(model = model_str, data = df_mg,
               missing = "fiml", estimator = "MLR")
  if (has_cluster) args$cluster <- "id_fam"
  if (multi_group) {
    args$group       <- "group"
    args$group.equal <- c("regressions", "lv.variances", "lv.covariances")
  }
  do.call(sem, args)
}

df_unrelated <- read.csv(paste0(sourceFileStem, "GN_scaled_selectunpaired.csv"))

subgroups <- list(
  list(label = "Whole_Sample",  data = df_sem,       include_sex = TRUE,  trios_only = FALSE),
  list(label = "All_Trios",     data = df_sem,       include_sex = TRUE,  trios_only = TRUE ),
  list(label = "All_Unrelated", data = df_unrelated, include_sex = TRUE,  trios_only = FALSE),
  list(label = "All_Females",   data = dplyr::filter(df_sem, sex1 == 0), include_sex = FALSE, trios_only = FALSE),
  list(label = "All_Males",     data = dplyr::filter(df_sem, sex1 == 1), include_sex = FALSE, trios_only = FALSE)
)

for (sg in subgroups) {
  sg_label <- sg$label
  sg_data  <- sg$data
  cat(sprintf("\n========== Subgroup: %s (N rows = %d) ==========\n", sg_label, nrow(sg_data)))

  for (i in seq_along(VarPairs)) {
    vp          <- VarPairs[[i]]
    trait_label <- vp[6]
    y_var       <- paste0(vp[1], "1")
    child_pgs   <- vp[2]
    mom_pgs     <- paste0(vp[2], "_fp")
    dad_pgs     <- paste0(vp[2], "_mp")

    covs <- if (sg$include_sex) c(vp[3], "sex1", "chiptype") else c(vp[3], "chiptype")

    missing_covs <- setdiff(covs, names(sg_data))
    if (length(missing_covs))
      stop(sprintf("[%s | %s] missing covariates: %s", sg_label, trait_label,
                   paste(missing_covs, collapse = ", ")))
    covs_str <- paste(covs, collapse = " + ")

    df_mg <- build_df_mg(sg_data, y_var, child_pgs, mom_pgs, dad_pgs, covs)
    if (sg$trios_only) df_mg <- df_mg %>% filter(group == "trios")

    n_var <- length(covs) + 4
    df_mg <- df_mg %>% group_by(group) %>% filter(n() >= n_var) %>% ungroup()

    cat(sprintf("  [%s | %d/%d] %s | N total = %d | trios = %d | dad_missing = %d | mom_missing = %d\n",
                sg_label, i, length(VarPairs), trait_label, nrow(df_mg),
                sum(df_mg$group=="trios"), sum(df_mg$group=="dad_missing"),
                sum(df_mg$group=="mom_missing")))

    model_pop_str <- paste0('
      mom_lat =~ 1*', mom_pgs, '
      dad_lat =~ 1*', dad_pgs, '
      mom_lat ~~ vg*mom_lat
      dad_lat ~~ vg*dad_lat
      mom_lat ~~ cov_am*dad_lat
      ', child_pgs, ' ~ 0.5*mom_lat + 0.5*dad_lat
      ', child_pgs, ' ~~ res_c*', child_pgs, '
      ', y_var, ' ~ b_pop*', child_pgs, ' + ', covs_str, '
      ', y_var, ' ~~ ve*', y_var, '
    ')
    model_full_str <- paste0('
      mom_lat =~ 1*', mom_pgs, '
      dad_lat =~ 1*', dad_pgs, '
      mom_lat ~~ vg*mom_lat
      dad_lat ~~ vg*dad_lat
      mom_lat ~~ cov_am*dad_lat
      ', child_pgs, ' ~ 0.5*mom_lat + 0.5*dad_lat
      ', child_pgs, ' ~~ res_c*', child_pgs, '
      ', y_var, ' ~ b_child*', child_pgs, ' + b_mom*mom_lat + b_dad*dad_lat + ', covs_str, '
      ', y_var, ' ~~ ve*', y_var, '
    ')

    fit_pop  <- fit_one(model_pop_str,  df_mg)
    fit_full <- fit_one(model_full_str, df_mg)

    if (!lavInspect(fit_pop,  "converged")) warning(sprintf("[%s | %s] Population model did not converge", sg_label, trait_label))
    if (!lavInspect(fit_full, "converged")) warning(sprintf("[%s | %s] Genetic nurture model did not converge",  sg_label, trait_label))

    res_pop  <- extract_param(fit_pop,  "b_pop",                        trait_label, "Pop")
    res_full <- extract_param(fit_full, c("b_child", "b_mom", "b_dad"), trait_label, "GN")

    res_combined <- bind_rows(res_pop, res_full)
    res_combined$Subgroup    <- sg_label
    res_combined$domain      <- vp[5]
    res_combined$n_trios     <- sum(df_mg$group == "trios")
    res_combined$n_total     <- nrow(df_mg)

    trio_rows <- df_mg$group == "trios"
    s_mom <- sd(df_mg[[mom_pgs]][trio_rows], na.rm = TRUE)
    s_dad <- sd(df_mg[[dad_pgs]][trio_rows], na.rm = TRUE)
    r_scale <- ifelse(res_combined$param == "b_mom", s_mom,
                      ifelse(res_combined$param == "b_dad", s_dad, 1))
    res_combined$est_r      <- res_combined$est      * r_scale
    res_combined$se_r       <- res_combined$se       * r_scale
    res_combined$ci_lower_r <- res_combined$ci_lower * r_scale
    res_combined$ci_upper_r <- res_combined$ci_upper * r_scale

    key <- paste(sg_label, i, sep = "__")
    all_results[[key]] <- res_combined

    fit_row <- bind_rows(extract_fit(fit_pop,  trait_label, "Pop"),
                         extract_fit(fit_full, trait_label, "GN"))
    fit_row$Subgroup <- sg_label
    fit_row$domain   <- vp[5]
    fit_results[[key]] <- fit_row

    cat("    Done.\n")
  }
}

df_results <- bind_rows(all_results)
df_fit     <- bind_rows(fit_results)

df_results <- df_results %>%
  group_by(param, model, domain, Subgroup) %>%
  mutate(pvalue_fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

df_fit_out <- df_fit %>%
  dplyr::rename(Unique_Label = trait, Domain = domain) %>%
  dplyr::relocate(Domain, Unique_Label, Subgroup)

make_block <- function(m, p, stem, is_parent) {
  d <- df_results[df_results$model == m & df_results$param == p, , drop = FALSE]
  out <- data.frame(trait = d$trait, Subgroup = d$Subgroup, stringsAsFactors = FALSE)
  out[[paste0(stem, "_Est")]]     <- d$est
  out[[paste0(stem, "_SE")]]      <- d$se
  out[[paste0(stem, "_LowerCI")]] <- d$ci_lower
  out[[paste0(stem, "_UpperCI")]] <- d$ci_upper
  out[[paste0(stem, "_P")]]       <- d$pvalue
  out[[paste0(stem, "_P_FDR")]]   <- d$pvalue_fdr
  if (is_parent) {
    out[[paste0(stem, "_r_Est")]]     <- d$est_r
    out[[paste0(stem, "_r_SE")]]      <- d$se_r
    out[[paste0(stem, "_r_LowerCI")]] <- d$ci_lower_r
    out[[paste0(stem, "_r_UpperCI")]] <- d$ci_upper_r
  }
  out
}

wide_spec <- list(
  list("Pop", "b_pop",   "Beta_Child_Pop",       FALSE),
  list("GN",  "b_child", "Beta_Child_OverallGN", FALSE),
  list("GN",  "b_mom",   "Beta_Mom_OverallGN",   TRUE ),
  list("GN",  "b_dad",   "Beta_Dad_OverallGN",   TRUE )
)

df_results_out <- df_results %>%
  distinct(trait, Subgroup, domain, n_trios, n_total)
for (sp in wide_spec)
  df_results_out <- left_join(df_results_out,
                              make_block(sp[[1]], sp[[2]], sp[[3]], sp[[4]]),
                              by = c("trait", "Subgroup"))

df_results_out <- df_results_out %>%
  rename(Domain = domain, Unique_Label = trait, N_used = n_total, N_trio = n_trios) %>%
  relocate(Domain, Unique_Label, Subgroup, N_used, N_trio)

wb <- createWorkbook()
addWorksheet(wb, "SEM")
writeData(wb, "SEM", df_results_out)
addWorksheet(wb, "SEM_ModelFit")
writeData(wb, "SEM_ModelFit", df_fit_out)
saveWorkbook(
  wb,
  file = paste0(outFileStem, "GN_SEM.xlsx"),
  overwrite = TRUE
)

write.csv(df_results_out,
          paste0(outFileStem, "GN_SEM.csv"),
          row.names = FALSE)
write.csv(df_fit_out,
          paste0(outFileStem, "GN_SEM_ModelFit.csv"),
          row.names = FALSE)

trait_levels_ordered <- rev(sapply(VarPairs, function(vp) vp[6]))
clean_trait_name <- function(x) {
  x <- gsub("^childhood\\s*|^adolescence\\s*|^adulthood\\s*", "", x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
trait_levels_cleaned <- unique(sapply(trait_levels_ordered, clean_trait_name))

df_plot <- df_results %>%
  filter(Subgroup == "Whole_Sample",
         param %in% c("b_mom", "b_dad", "b_pop", "b_child")) %>%
  mutate(
    param_label = case_when(
      param == "b_mom"   ~ "Maternal PGS (Genetic Nurture)",
      param == "b_dad"   ~ "Paternal PGS (Genetic Nurture)",
      param == "b_pop"   ~ "Child PGS (Population)",
      param == "b_child" ~ "Child PGS (Direct)",
      TRUE ~ param
    ),
    param_label = factor(param_label,
                         levels = c("Paternal PGS (Genetic Nurture)", "Maternal PGS (Genetic Nurture)",
                                    "Child PGS (Direct)", "Child PGS (Population)")),

    trait_clean = gsub("^childhood\\s*|^adolescence\\s*|^adulthood\\s*", "", trait),
    trait_clean = paste0(toupper(substr(trait_clean, 1, 1)), substr(trait_clean, 2, nchar(trait_clean))),
    trait_ordered = factor(trait_clean, levels = trait_levels_cleaned),

    star_label = fdr_stars(pvalue_fdr)
  )

color_map <- SERIES4_COL

p_forest <- ggplot(df_plot, aes(x = est_r, y = trait_ordered, color = param_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower_r, xmax = ci_upper_r), orientation = "y",
                width = 0.25, linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
  geom_text(
    aes(x = ci_upper_r, label = star_label),
    hjust = -0.3, vjust = 0.75, size = 4.5, show.legend = FALSE,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~ domain, ncol = 2, scales = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  scale_color_manual(values = color_map, name = NULL) +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_x_continuous(name = X_BETA, breaks = BRK_BETA) +
  labs(
    title    = "Child PGS effects and parental genetic nurture effects: three-group SEM",
  ) +
  theme_gn(grid = "x") +
  theme(axis.title.y  = element_blank(),
        panel.spacing = unit(1.5, "lines"))

save_gn(p_forest, "SEM_forest_GN_vs_Pop_2x2", width = 18, height = 18)

FIG_MAT <- MAT_COL
FIG_PAT <- PAT_COL
FIG_SEX <- unname(SEX_COL)
names(FIG_SEX) <- c("Female", "Male")

estimand_map <- data.frame(
  model    = c("Pop",    "GN",      "GN",    "GN"),
  param    = c("b_pop",  "b_child", "b_mom", "b_dad"),
  estimand = c("Child Population", "Child | Both Parents",
               "Maternal Genetic Nurture", "Paternal Genetic Nurture"),
  stringsAsFactors = FALSE
)

GN_FOREST_XLIM <- {
  d <- df_results %>%
    filter(model == "GN", param %in% c("b_mom", "b_dad"),
           Subgroup %in% c("All_Females", "All_Males",
                           "All_Trios", "All_Unrelated", "Whole_Sample"))
  lo <- min(d$ci_lower_r, na.rm = TRUE); hi <- max(d$ci_upper_r, na.rm = TRUE)
  pad <- (hi - lo) * 0.08
  c(lo - pad, hi + pad)
}

df_sex <- df_results %>%
  filter(Subgroup %in% c("All_Females", "All_Males")) %>%
  inner_join(estimand_map, by = c("model", "param")) %>%
  select(trait, domain, estimand, Subgroup, est, se) %>%
  pivot_wider(names_from = Subgroup, values_from = c(est, se)) %>%
  mutate(
    z = (est_All_Females - est_All_Males) /
      sqrt(se_All_Females^2 + se_All_Males^2),
    p_sex = 2 * pnorm(-abs(z))
  ) %>%
  group_by(estimand, domain) %>%
  mutate(p_sex_fdr = p.adjust(p_sex, method = "fdr")) %>%
  ungroup() %>%
  mutate(estimand = factor(estimand, levels = estimand_map$estimand),
         sig_sex  = p_sex_fdr < 0.05)

if (anyNA(df_sex$est_All_Females) || anyNA(df_sex$est_All_Males))
  stop("Missing All_Females / All_Males estimates; check the converged column in GN_SEM_ModelFit.csv")

sex_ann <- df_sex %>%
  group_by(estimand) %>%
  summarise(r      = cor(est_All_Females, est_All_Males),
            mean_d = mean(est_All_Females - est_All_Males),
            n_sig  = sum(sig_sex), n = n(), .groups = "drop") %>%
  mutate(lab = sprintf("r = %.2f\nmean \u0394 (F\u2212M) = %+.3f\nsex diff, FDR q<.05 : %d/%d",
                       r, mean_d, n_sig, n))

plot_sex_forest <- function(param_lbl, estimand_lbl, file_stem) {
  d <- df_results %>%
    filter(Subgroup %in% c("All_Females", "All_Males"), model == "GN", param == param_lbl) %>%
    mutate(
      trait_clean   = clean_trait_name(trait),
      trait_ordered = factor(trait_clean, levels = trait_levels_cleaned),
      Sex   = factor(ifelse(Subgroup == "All_Females", "Female", "Male"), levels = c("Male", "Female")),
      star  = fdr_stars(pvalue_fdr)
    )

  dg <- position_dodge(width = 0.5)

  p <- ggplot(d, aes(x = est_r, y = trait_ordered, colour = Sex)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lower_r, xmax = ci_upper_r), orientation = "y",
                  width = 0, linewidth = 0.75, position = dg) +
    geom_point(size = 2.0, position = dg) +
    geom_text(aes(x = ci_upper_r, label = star), hjust = -0.35, vjust = 0.75,
              size = 5, show.legend = FALSE, position = dg) +
    facet_wrap(~ domain, ncol = 2, scales = "free_y",
               labeller = label_wrap_gen(width = 22)) +
    scale_colour_manual(values = FIG_SEX, name = NULL,
                        guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(name = X_BETA, breaks = BRK_BETA, limits = GN_FOREST_XLIM) +
    labs(title = paste0(estimand_lbl, " in female vs. male offspring: three-group SEM")) +
    theme_gn(grid = "x") +
    theme(axis.title.y  = element_blank(),
          panel.spacing = unit(1.5, "lines"))

  save_gn(p, file_stem, width = 16, height = 14)
}

plot_sex_forest("b_mom", "Maternal genetic nurture", "SEM_sex_forest_maternal")
plot_sex_forest("b_dad", "Paternal genetic nurture", "SEM_sex_forest_paternal")

SAMPLE_SGS <- c("All_Trios", "All_Unrelated", "Whole_Sample")

plot_sample_forest <- function(param_lbl, estimand_lbl, hi_col, file_stem) {
  d <- df_results %>%
    filter(Subgroup %in% SAMPLE_SGS, model == "GN", param == param_lbl) %>%
    mutate(
      trait_clean   = clean_trait_name(trait),
      trait_ordered = factor(trait_clean, levels = trait_levels_cleaned),
      Subgroup      = factor(Subgroup, levels = rev(SAMPLE_SGS)),
      star          = fdr_stars(pvalue_fdr)
    )

  pal <- sample_ramp(hi_col)
  dg  <- position_dodge(width = 0.68)

  p <- ggplot(d, aes(x = est_r, y = trait_ordered, colour = Subgroup)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lower_r, xmax = ci_upper_r), orientation = "y",
                  width = 0, linewidth = 0.75, position = dg) +
    geom_point(size = 2.0, position = dg) +
    geom_text(aes(x = ci_upper_r, label = star), hjust = -0.35, vjust = 0.75,
              size = 5, show.legend = FALSE, position = dg) +
    facet_wrap(~ domain, ncol = 2, scales = "free_y",
               labeller = label_wrap_gen(width = 22)) +
    scale_colour_manual(values = pal, name = NULL,
                        labels = c(All_Trios     = "Complete Trios",
                                   All_Unrelated = "Three-Group SEM (Unrelated)",
                                   Whole_Sample  = "Three-Group SEM (Whole Sample)"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_x_continuous(name = X_BETA, breaks = BRK_BETA, limits = GN_FOREST_XLIM) +
    labs(title = paste0(estimand_lbl, ": sample-composition sensitivity")) +
    theme_gn(grid = "x") +
    theme(axis.title.y  = element_blank(),
          panel.spacing = unit(1.5, "lines"))

  save_gn(p, file_stem, width = 16, height = 14)
}

plot_sample_forest("b_mom", "Maternal genetic nurture", FIG_MAT, "SEM_sample_forest_maternal")
plot_sample_forest("b_dad", "Paternal genetic nurture", FIG_PAT, "SEM_sample_forest_paternal")

cat("\n--- Sex differences in genetic nurture (z test on SEM coefficients) ---\n")
print(as.data.frame(sex_ann[, c("estimand", "r", "mean_d", "n_sig", "n")]), row.names = FALSE)
if (any(df_sex$sig_sex)) {
  cat("\nTraits with a sex difference surviving FDR:\n")
  print(as.data.frame(df_sex %>% filter(sig_sex) %>%
                        select(estimand, domain, trait, est_All_Females, est_All_Males, p_sex_fdr)),
        row.names = FALSE)
}

cat("\n--- Median SE ratio and number of FDR-significant estimates relative to complete trios ---\n")
print(as.data.frame(
  df_results %>%
    filter(Subgroup %in% SAMPLE_SGS) %>%
    inner_join(estimand_map, by = c("model", "param")) %>%
    group_by(estimand, Subgroup) %>%
    summarise(median_SE = median(se), n_sig_FDR = sum(pvalue_fdr < 0.05), .groups = "drop") %>%
    group_by(estimand) %>%
    mutate(SE_ratio_vs_Trios = round(median_SE / median_SE[Subgroup == "All_Trios"], 3),
           median_SE = round(median_SE, 4)) %>%
    ungroup() %>%
    arrange(match(estimand, estimand_map$estimand), match(Subgroup, SAMPLE_SGS))
), row.names = FALSE)

reg_all <- read.csv(paste0(regFileStem, "GN_Reg.csv"), stringsAsFactors = FALSE) %>%
  filter(Subgroup == "All_Trios")
sem_all <- df_results %>% filter(Subgroup == "All_Trios")

compare_spec <- data.frame(
  model    = c("Pop",           "GN",           "GN",                 "GN"),
  param    = c("b_pop",         "b_child",      "b_mom",              "b_dad"),
  stem     = c("Beta_Child_Pop", "Beta_Child_OverallGN",
               "Beta_Mom_OverallGN", "Beta_Dad_OverallGN"),
  panel    = c("Child PGS (Population)", "Child PGS (Direct)",
               "Maternal PGS (Genetic Nurture)", "Paternal PGS (Genetic Nurture)"),
  stringsAsFactors = FALSE
)

cmp <- bind_rows(lapply(seq_len(nrow(compare_spec)), function(i) {
  sp <- compare_spec[i, ]
  s  <- sem_all %>% filter(model == sp$model, param == sp$param) %>%
    select(trait, domain, sem_est = est, sem_se = se,
           sem_lo = ci_lower, sem_hi = ci_upper)
  r  <- reg_all %>%
    select(trait = Unique_Label,
           reg_est = all_of(paste0(sp$stem, "_Est")),
           reg_se  = all_of(paste0(sp$stem, "_SE")),
           reg_lo  = all_of(paste0(sp$stem, "_LowerCI")),
           reg_hi  = all_of(paste0(sp$stem, "_UpperCI")))
  inner_join(s, r, by = "trait") %>% mutate(panel = sp$panel)
})) %>%
  mutate(panel = factor(panel, levels = compare_spec$panel),
         delta = sem_est - reg_est,
         se_ratio = sem_se / reg_se)

n_expect <- length(VarPairs) * nrow(compare_spec)
if (nrow(cmp) != n_expect)
  stop(sprintf("SEM and regression results do not match: %d rows after the join, expected %d (%d traits x %d estimands)",
               nrow(cmp), n_expect, length(VarPairs), nrow(compare_spec)))
if (anyNA(cmp$sem_est) || anyNA(cmp$reg_est) || anyNA(cmp$sem_lo) || anyNA(cmp$reg_lo))
  stop("Missing estimates or CIs; check the converged column in GN_SEM_ModelFit.csv")

p_semvsreg <- ggplot(cmp, aes(x = reg_est, y = sem_est, colour = domain)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(xmin = reg_lo, xmax = reg_hi),
                orientation = "y", width = 0, linewidth = 0.35, alpha = 0.25) +
  geom_errorbar(aes(ymin = sem_lo, ymax = sem_hi),
                orientation = "x", width = 0, linewidth = 0.35, alpha = 0.25) +
  geom_point(size = 2.3, alpha = 0.85) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = DOM_COL, labels = DOM_SHORT, name = NULL) +
  scale_x_continuous(breaks = BRK_BETA) +
  scale_y_continuous(breaks = BRK_BETA) +
  coord_equal() +
  labs(x = "Regression: Beta (95% CI)",
       y = "Three-group SEM: Beta (95% CI)",
       title = "Three-group SEM and regression give the same estimates in complete trios",
  ) +
  theme_gn(grid = "both")

save_gn(p_semvsreg, "SEM_vs_Reg_CompleteTrios", width = 20, height = 6.4)

cat("\n--- Three-group SEM versus regression in complete trios ---\n")
print(as.data.frame(
  cmp %>% group_by(panel) %>%
    summarise(n              = n(),
              max_abs_delta  = max(abs(delta)),
              max_rel_delta  = max(abs(delta) / pmax(abs(reg_est), 1e-12)),
              median_SE_ratio = median(se_ratio),
              max_SE_ratio    = max(se_ratio),
              .groups = "drop")
), row.names = FALSE)
