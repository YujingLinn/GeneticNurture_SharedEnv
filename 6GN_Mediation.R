library(lavaan)
library(dplyr)
library(tibble)
library(openxlsx)
library(ggplot2)
library(tidyr)
library(stringr)

source("0GN_Theme.R")

sourceFileStem <- "GN_Data/"
outFileStem    <- "GN_Results/Mediation/"
dir.create(outFileStem, recursive = TRUE, showWarnings = FALSE)

source("0GN_VarPairs.R")

dat_scaled <- read.csv(paste0(sourceFileStem, "GN_scaled.csv"))

ses_vars_all <- c("ases", "gses", "pses", "u1pses")
missing_ses  <- setdiff(ses_vars_all, names(dat_scaled))
if (length(missing_ses) > 0) {
  stop(sprintf("SES columns missing from GN_scaled.csv: %s -- check 0GN_Clean.R output before proceeding.",
               paste(missing_ses, collapse = ", ")))
}

ses_wave_matched <- function(age_num) {
  if (age_num <= 14) return("gses")
  if (age_num <= 18) return("pses")
  return("u1pses")
}
ses_wave_fixed <- function(age_num) "ases"

extract_age_num <- function(unique_label) as.numeric(gsub(".*age ([0-9]+).*", "\\1", unique_label))

run_mediation <- function(data, var_pairs_list, ses_wave_fun, scheme_label) {

  out_rows <- list()

  for (vp in var_pairs_list) {
    outcome_var <- paste0(vp[1], "1")
    pgs_base    <- vp[2]
    age_var     <- vp[3]
    domain      <- vp[5]
    unique_lbl  <- vp[6]

    child_pgs  <- pgs_base
    mother_pgs <- paste0(pgs_base, "_fp")
    father_pgs <- paste0(pgs_base, "_mp")
    ses_var    <- ses_wave_fun(extract_age_num(unique_lbl))

    vars_needed <- c(outcome_var, child_pgs, mother_pgs, father_pgs, ses_var,
                     age_var, "sex1", "chiptype", "id_fam")
    missing_cols <- setdiff(vars_needed, names(data))
    if (length(missing_cols) > 0) {
      cat(sprintf("[%s | %s] SKIPPING: missing columns %s\n", scheme_label, unique_lbl,
                  paste(missing_cols, collapse = ", ")))
      next
    }

    df_clean <- data[complete.cases(data[, vars_needed]), vars_needed]
    if (nrow(df_clean) < 100) {
      cat(sprintf("[%s | %s] SKIPPING: N = %d too small\n", scheme_label, unique_lbl, nrow(df_clean)))
      next
    }

    covar_rhs <- paste(age_var, "sex1", "as.factor(chiptype)", sep = " + ")
    resid_of <- function(varname) {
      f <- as.formula(paste(varname, "~", covar_rhs))
      residuals(lm(f, data = df_clean))
    }
    df_clean$outcome_adjCovar     <- resid_of(outcome_var)
    df_clean$SES_adjCovar         <- resid_of(ses_var)
    df_clean$childPGS_adjCovar    <- resid_of(child_pgs)
    df_clean$maternalPGS_adjCovar <- resid_of(mother_pgs)
    df_clean$paternalPGS_adjCovar <- resid_of(father_pgs)

    model_str <- '
      outcome_adjCovar ~ c1*childPGS_adjCovar + c2*maternalPGS_adjCovar + c3*paternalPGS_adjCovar + b*SES_adjCovar
      SES_adjCovar     ~ a1*childPGS_adjCovar + a2*maternalPGS_adjCovar + a3*paternalPGS_adjCovar
      indirect_child := a1*b
      indirect_mom   := a2*b
      indirect_dad   := a3*b
      total_child    := c1 + indirect_child
      total_mom      := c2 + indirect_mom
      total_dad      := c3 + indirect_dad
    '

    fit <- sem(model_str, data = df_clean, cluster = "id_fam", estimator = "MLR")
    pe  <- parameterEstimates(fit)

    pick <- function(lbl) {
      row <- pe[pe$label == lbl, ]
      if (nrow(row) == 0) return(c(Est = NA_real_, SE = NA_real_, P = NA_real_))
      c(Est = row$est[1], SE = row$se[1], P = row$pvalue[1])
    }

    labels <- c("a1", "a2", "a3", "b", "c1", "c2", "c3",
                "indirect_child", "indirect_mom", "indirect_dad",
                "total_child", "total_mom", "total_dad")
    vals <- setNames(lapply(labels, pick), labels)

    s_child <- sd(df_clean$childPGS_adjCovar)
    s_mat   <- sd(df_clean$maternalPGS_adjCovar)
    s_pat   <- sd(df_clean$paternalPGS_adjCovar)
    s_ses   <- sd(df_clean$SES_adjCovar)
    r_scale <- c(a1 = s_child, a2 = s_mat, a3 = s_pat, b = s_ses,
                 c1 = s_child, c2 = s_mat, c3 = s_pat,
                 indirect_child = s_child, indirect_mom = s_mat, indirect_dad = s_pat,
                 total_child = s_child, total_mom = s_mat, total_dad = s_pat)

    row_out <- tibble(
      Domain = domain, Unique_Label = unique_lbl, Subgroup = scheme_label,
      SES_Var = ses_var, N_used = nrow(df_clean), N_fam = length(unique(df_clean$id_fam))
    )
    for (lbl in labels) {
      row_out[[paste0(lbl, "_Est")]]   <- vals[[lbl]]["Est"]
      row_out[[paste0(lbl, "_SE")]]    <- vals[[lbl]]["SE"]
      row_out[[paste0(lbl, "_P")]]     <- vals[[lbl]]["P"]
      row_out[[paste0(lbl, "_r_Est")]] <- vals[[lbl]]["Est"] * r_scale[[lbl]]
      row_out[[paste0(lbl, "_r_SE")]]  <- vals[[lbl]]["SE"]  * r_scale[[lbl]]
    }
    row_out$Atten_Pct_Mom <- 100 * vals$indirect_mom["Est"] / vals$total_mom["Est"]
    row_out$Atten_Pct_Dad <- 100 * vals$indirect_dad["Est"] / vals$total_dad["Est"]

    out_rows[[length(out_rows) + 1]] <- row_out
    cat(sprintf("[%s | %s] N=%d, SES=%s, done.\n", scheme_label, unique_lbl, nrow(df_clean), ses_var))
  }

  bind_rows(out_rows)
}

res_ases    <- run_mediation(dat_scaled, VarPairs, ses_wave_fixed,   "SES_Fixed_ases")
res_matched <- run_mediation(dat_scaled, VarPairs, ses_wave_matched, "SES_Matched")

SES_Mediation <- bind_rows(res_ases, res_matched)

path_labels <- c("a1", "a2", "a3", "b", "c1", "c2", "c3",
                 "indirect_child", "indirect_mom", "indirect_dad",
                 "total_child", "total_mom", "total_dad")
fdr_cols <- paste0(path_labels, "_P")
for (col in fdr_cols) {
  fdr_col <- paste0(col, "_FDR")
  SES_Mediation <- SES_Mediation %>%
    group_by(Subgroup, Domain) %>%
    mutate("{fdr_col}" := p.adjust(.data[[col]], method = "fdr")) %>%
    ungroup() %>%
    relocate(all_of(fdr_col), .after = all_of(col))
}

write.csv(SES_Mediation, file.path(outFileStem, "SES_Mediation.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "SES_Mediation")
writeData(wb, "SES_Mediation", SES_Mediation)
freezePane(wb, "SES_Mediation", firstRow = TRUE)
saveWorkbook(wb, file.path(outFileStem, "SES_Mediation.xlsx"), overwrite = TRUE)


GN_Reg_ref <- read.csv(file.path(outFileStem, "..", "Reg", "GN_Reg.csv")) %>%
  filter(Subgroup == "All_Trios") %>%
  select(Unique_Label,
         Beta_Mom_OverallGN_r_Est, Beta_Mom_OverallGN_r_SE, Beta_Mom_OverallGN_P_FDR,
         Beta_Dad_OverallGN_r_Est, Beta_Dad_OverallGN_r_SE, Beta_Dad_OverallGN_P_FDR)

trait_order  <- sapply(VarPairs, function(x) x[6])
domain_order <- c("Cognitive Abilities", "Education Achievement & Attainment",
                  "Height & BMI", "Psychopathology")

build_decomp_data <- function(parent = c("mom", "dad")) {
  parent <- match.arg(parent)
  if (parent == "mom") {
    dat <- SES_Mediation %>% filter(Subgroup == "SES_Fixed_ases") %>%
      select(Unique_Label, Domain,
             est_indirect = indirect_mom_r_Est, se_indirect = indirect_mom_r_SE,
             p_indirect   = indirect_mom_P_FDR,
             est_total    = total_mom_r_Est,    se_total    = total_mom_r_SE,
             p_total      = total_mom_P_FDR) %>%
      left_join(GN_Reg_ref %>% select(Unique_Label,
                                      est_reg = Beta_Mom_OverallGN_r_Est,
                                      se_reg  = Beta_Mom_OverallGN_r_SE,
                                      p_reg   = Beta_Mom_OverallGN_P_FDR),
                by = "Unique_Label")
    base_col <- MAT_COL
    lbl_reg  <- "Maternal PGS (Genetic Nurture, Regression Model)"
    lbl_tot  <- "Total Maternal Genetic Nurture (Mediation Model)"
    lbl_ind  <- "SES-Mediated Genetic Nurture"
    parent_label <- "Maternal"
  } else {
    dat <- SES_Mediation %>% filter(Subgroup == "SES_Fixed_ases") %>%
      select(Unique_Label, Domain,
             est_indirect = indirect_dad_r_Est, se_indirect = indirect_dad_r_SE,
             p_indirect   = indirect_dad_P_FDR,
             est_total    = total_dad_r_Est,    se_total    = total_dad_r_SE,
             p_total      = total_dad_P_FDR) %>%
      left_join(GN_Reg_ref %>% select(Unique_Label,
                                      est_reg = Beta_Dad_OverallGN_r_Est,
                                      se_reg  = Beta_Dad_OverallGN_r_SE,
                                      p_reg   = Beta_Dad_OverallGN_P_FDR),
                by = "Unique_Label")
    base_col <- PAT_COL
    lbl_reg  <- "Paternal PGS (Genetic Nurture, Regression Model)"
    lbl_tot  <- "Total Paternal Genetic Nurture (Mediation Model)"
    lbl_ind  <- "SES-Mediated Genetic Nurture"
    parent_label <- "Paternal"
  }

  pdat <- dat %>%
    pivot_longer(cols = -c(Unique_Label, Domain),
                 names_to = c(".value", "Model"),
                 names_pattern = "(est|se|p)_(.*)") %>%
    mutate(
      lo = est - 1.96 * se,
      hi = est + 1.96 * se,
      Model = dplyr::recode(Model, indirect = lbl_ind, total = lbl_tot, reg = lbl_reg),
      Model = factor(Model, levels = c(lbl_ind, lbl_tot, lbl_reg)),
      trait_clean = cap_first(Unique_Label),
      Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
      trait_clean  = factor(trait_clean,
                            levels = unique(trait_clean[order(match(Unique_Label, levels(Unique_Label)))])),
      Domain = factor(Domain, levels = domain_order),
      Stars  = fdr_stars(p)
    )

  list(pdat = pdat, base_col = base_col, lbl_reg = lbl_reg, lbl_tot = lbl_tot,
       lbl_ind = lbl_ind, parent_label = parent_label)
}

mom_prep <- build_decomp_data("mom")
dad_prep <- build_decomp_data("dad")

all_bounds  <- c(mom_prep$pdat$lo, mom_prep$pdat$hi, dad_prep$pdat$lo, dad_prep$pdat$hi)
bounds_rng  <- range(all_bounds, na.rm = TRUE)
pad         <- diff(bounds_rng) * 0.10
shared_xlim <- bounds_rng + c(-pad, pad)

render_decomp_plot <- function(prep, xlim) {
  pdat <- prep$pdat; base_col <- prep$base_col
  lbl_reg <- prep$lbl_reg; lbl_tot <- prep$lbl_tot; lbl_ind <- prep$lbl_ind
  parent_label <- prep$parent_label

  model_colors <- setNames(c(lighten_gn(base_col), base_col, OVR_COL), c(lbl_ind, lbl_tot, lbl_reg))
  dodge_width <- 0.65

  ggplot(pdat, aes(x = est, y = trait_clean, color = Model)) +
    geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey50") +
    geom_linerange(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = dodge_width), linewidth = 0.8) +
    geom_point(size = 3, position = position_dodge(width = dodge_width)) +
    geom_text(aes(x = hi + diff(xlim) * 0.02, label = Stars),
              position = position_dodge(width = dodge_width),
              vjust = 0.75, size = 4.5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = model_colors, breaks = c(lbl_reg, lbl_tot, lbl_ind)) +
    scale_x_continuous(breaks = BRK_BETA, limits = xlim) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    facet_wrap(~ Domain, scales = "free_y", ncol = 2) +
    labs(
      title = paste0(parent_label, " genetic nurture effect: extent mediated by family SES"),
      x     = X_BETA,
      y     = NULL,
      color = NULL
    ) +
    theme_gn(grid = "x") +
    guides(color = guide_legend(nrow = 1))
}

plot_mom_decomp <- render_decomp_plot(mom_prep, shared_xlim)
plot_dad_decomp <- render_decomp_plot(dad_prep, shared_xlim)

calc_height <- 4.0 + ((length(unique(SES_Mediation$Unique_Label)) / 2) * 0.8)
save_gn(plot_mom_decomp, "SES_Mediation_Decomposition_Mom", width = 24, height = calc_height)
save_gn(plot_dad_decomp, "SES_Mediation_Decomposition_Dad", width = 24, height = calc_height)
