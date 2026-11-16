sourceFileStem <- "GN_Data"
outFileStem    <- file.path("GN_Results", "Data_Checking")
dir.create(outFileStem, recursive = TRUE, showWarnings = FALSE)

library(dplyr); library(tidyr); library(purrr); library(tibble)
library(moments); library(ggplot2); library(openxlsx)

df_full <- read.csv(file.path(sourceFileStem, "GN_raw.csv"))
df      <- read.csv(file.path(sourceFileStem, "GN_raw_selectunpaired.csv"))
source("0GN_VarPairs.R")
source("0GN_Theme.R")

pheno_vars   <- paste0(sapply(VarPairs, function(x) x[1]), "1")
pheno_labels <- sapply(VarPairs, function(x) x[6])
pheno_types  <- sapply(VarPairs, function(x) x[7])

ses_vars   <- c("ases", "gses", "pses", "u1pses")
ses_labels <- c("SES (Birth)", "SES (Age 7)", "SES (Age 16)", "SES (Age 21)")

child_pgs_vars  <- unique(sapply(VarPairs, function(x) x[2]))
mother_pgs_vars <- paste0(child_pgs_vars, "_fp")
father_pgs_vars <- paste0(child_pgs_vars, "_mp")
pgs_short <- c(
  "gFactor_sumstats_noMAF_IaddedN_LDpred2.txt_pred_auto"        = "g Factor",
  "EA4_no23andme_effN_allsnp_Okbay2022_ldpred2_pred_auto"       = "Edu Attain",
  "Height_Yengo2022_ldpred2_pred_auto"                          = "Height",
  "childhoodBMI_Vogelezang2020_YLprocessed.txt_pred_auto"       = "Childhood BMI",
  "BMI_Giant2018_pred_auto"                                     = "BMI",
  "PTSD_Nievergelt2024_EURv4_cont_ldpred2_noChrX.txt_pred_auto" = "PTSD",
  "ANX_2026_Skelton_Cont_PGC.txt_pred_auto"                     = "Anxiety",
  "MDD_Adams2025_no23andMe_EUR_casecontrol_ldpred2_pred_auto"   = "MDD",
  "AUDIT_UKB_2018_AJP_with_bp_filtered_LDpred2_pred_auto"       = "Alcohol",
  "Anorexia_2026_Termorshuizen_CaseControl_PGC_LDPred2.txt_pred_auto" = "Anorexia",
  "ADHD2022_iPSYCH_deCODE_PGC_meta_pred_auto"                   = "ADHD",
  "Bipolar_OConnell2025_no23andMe_EUR_cont_ldpred2_pred_auto"   = "Bipolar",
  "ASD_iPSYCH_PGC_2017_pred_auto"                               = "Autism"
)
child_pgs_labels  <- paste0("PGS child: ",  pgs_short[child_pgs_vars])
mother_pgs_labels <- paste0("PGS mother: ", pgs_short[child_pgs_vars])
father_pgs_labels <- paste0("PGS father: ", pgs_short[child_pgs_vars])

var_label_map <- c(setNames(pheno_labels, pheno_vars),
                   setNames(ses_labels,   ses_vars),
                   setNames(child_pgs_labels,  child_pgs_vars),
                   setNames(mother_pgs_labels, mother_pgs_vars),
                   setNames(father_pgs_labels, father_pgs_vars))
var_type_map  <- c(setNames(pheno_types, pheno_vars),
                   setNames(rep(NA_character_, length(ses_vars)),         ses_vars),
                   setNames(rep("Child PGS",  length(child_pgs_vars)),    child_pgs_vars),
                   setNames(rep("Mother PGS", length(mother_pgs_vars)),   mother_pgs_vars),
                   setNames(rep("Father PGS", length(father_pgs_vars)),   father_pgs_vars))

pgs_corr_map <- do.call(rbind, lapply(child_pgs_vars, function(v) data.frame(
  Var   = c(v, paste0(v, "_fp"), paste0(v, "_mp"),
            paste0(v, "_noPC"), paste0(v, "_fp_noPC"), paste0(v, "_mp_noPC")),
  Label = c(paste0(c("PGS child (.PC): ",   "PGS mother (.PC): ",   "PGS father (.PC): "),   pgs_short[v]),
            paste0(c("PGS child (noPC): ", "PGS mother (noPC): ", "PGS father (noPC): "), pgs_short[v])),
  stringsAsFactors = FALSE)))

corr_map <- bind_rows(
  pgs_corr_map,
  data.frame(Var = pheno_vars,        Label = paste0(pheno_types, " | ", pheno_labels), stringsAsFactors = FALSE),
  data.frame(
    Var   = c(paste0("PC", 1:10),            paste0("PC", 1:10, "_fp"),        paste0("PC", 1:10, "_mp")),
    Label = c(paste0("PC", 1:10, " (child)"), paste0("PC", 1:10, " (mother)"), paste0("PC", 1:10, " (father)")),
    stringsAsFactors = FALSE),
  data.frame(Var = ses_vars,          Label = ses_labels,                               stringsAsFactors = FALSE))
corr_map <- corr_map[!duplicated(corr_map$Var), ]

corr_present <- corr_map$Var %in% names(df)
corr_map <- corr_map[corr_present, ]

corr_matrix <- cor(df[, corr_map$Var], use = "pairwise.complete.obs")
rownames(corr_matrix) <- corr_map$Label; colnames(corr_matrix) <- corr_map$Label

write.csv(corr_matrix, file.path(outFileStem, "GN_CorrTable.csv"), row.names = TRUE)

pair_cor <- function(data, va, vb) {
  if (!all(c(va, vb) %in% names(data))) return(c(r = NA_real_, N = 0, lo = NA_real_, hi = NA_real_))
  a <- data[[va]]; b <- data[[vb]]; ok <- !is.na(a) & !is.na(b); n <- sum(ok)
  if (n < 4) return(c(r = if (n >= 3) cor(a[ok], b[ok]) else NA_real_, N = n, lo = NA_real_, hi = NA_real_))
  ct <- cor.test(a[ok], b[ok])
  c(r = unname(ct$estimate), N = n, lo = ct$conf.int[1], hi = ct$conf.int[2])
}
family_cor <- function(data) {
  map_dfr(seq_along(child_pgs_vars), function(k) {
    cm <- pair_cor(data, child_pgs_vars[k],  mother_pgs_vars[k])
    cf <- pair_cor(data, child_pgs_vars[k],  father_pgs_vars[k])
    mf <- pair_cor(data, mother_pgs_vars[k], father_pgs_vars[k])
    tibble(Trait = unname(pgs_short[child_pgs_vars[k]]),
           r_MotherChild = unname(cm["r"]), r_MotherChild_lowerCI = unname(cm["lo"]), r_MotherChild_upperCI = unname(cm["hi"]),
           r_FatherChild = unname(cf["r"]), r_FatherChild_lowerCI = unname(cf["lo"]), r_FatherChild_upperCI = unname(cf["hi"]),
           r_MotherFather_AM = unname(mf["r"]), r_MotherFather_AM_lowerCI = unname(mf["lo"]), r_MotherFather_AM_upperCI = unname(mf["hi"]))
  })
}

pgs_disp <- c("g Factor" = "Cognitive ability (g)", "Edu Attain" = "Years of Schooling",
              "Height" = "Height", "Childhood BMI" = "Childhood BMI", "BMI" = "Adult BMI",
              "PTSD" = "PTSD", "Anxiety" = "Anxiety", "MDD" = "Major Depression",
              "Alcohol" = "Alcohol Use (AUDIT)", "Anorexia" = "Eating Disorder (Anorexia)",
              "ADHD" = "ADHD", "Bipolar" = "Bipolar Disorder", "Autism" = "Autism")
fmt_ci <- function(r, lo, hi) ifelse(is.na(r), NA_character_, sprintf("%.2f (%.2f, %.2f)", r, lo, hi))
fmt_family_cor <- function(fc) {
  fc %>% filter(Trait != "N (pairs)") %>% transmute(
    Trait = unname(ifelse(Trait %in% names(pgs_disp), pgs_disp[Trait], Trait)),
    `r_Mother-Child (95% CI)`  = fmt_ci(r_MotherChild,     r_MotherChild_lowerCI,     r_MotherChild_upperCI),
    `r_Father-Child (95% CI)`  = fmt_ci(r_FatherChild,     r_FatherChild_lowerCI,     r_FatherChild_upperCI),
    `r_Mother-Father (95% CI)` = fmt_ci(r_MotherFather_AM, r_MotherFather_AM_lowerCI, r_MotherFather_AM_upperCI))
}

calc_desc <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(tibble(N = 0, Mean = NA_real_, SD = NA_real_, Min = NA_real_,
                  Max = NA_real_, Skewness = NA_real_, Kurtosis = NA_real_))
  tibble(N = length(x), Mean = mean(x), SD = sd(x), Min = min(x), Max = max(x),
         Skewness = if (length(x) > 2) skewness(x) else NA_real_,
         Kurtosis = if (length(x) > 3) kurtosis(x) else NA_real_)
}
age_of <- setNames(sapply(VarPairs, function(x) x[3]), pheno_vars)

desc_block <- function(data, vars, types, labels, category) {
  pmap_dfr(list(vars, types, labels), function(v, ty, lb) {
    if (!v %in% names(data)) return(NULL)
    sub <- data[!is.na(data[[v]]), ]
    age_var <- unname(age_of[v]); am <- NA_real_; asd <- NA_real_
    if (!is.na(age_var) && age_var %in% names(sub)) {
      av <- sub[[age_var]][!is.na(sub[[age_var]])]
      if (length(av) > 0) am  <- mean(av)
      if (length(av) > 1) asd <- sd(av)
    }
    pf  <- if ("sex1"  %in% names(sub)) mean(sub$sex1  == 0, na.rm = TRUE) else NA_real_
    pdz <- if ("zygos" %in% names(sub)) mean(sub$zygos == 2, na.rm = TRUE) else NA_real_
    bind_cols(tibble(Category = category, Trait_Type = ty, Label = lb, Variable = v,
                     Prop_Female = pf, Prop_DZ = pdz, Age_Mean = am, Age_SD = asd),
              calc_desc(data[[v]]))
  })
}
desc_sample <- function(data) {
  bind_rows(
    desc_block(data, pheno_vars,      pheno_types,                                 pheno_labels,      "Phenotype"),
    desc_block(data, ses_vars,        rep(NA_character_, length(ses_vars)),        ses_labels,        "SES"),
    desc_block(data, child_pgs_vars,  rep(NA_character_, length(child_pgs_vars)),  child_pgs_labels,  "Child PGS"),
    desc_block(data, mother_pgs_vars, rep(NA_character_, length(mother_pgs_vars)), mother_pgs_labels, "Mother PGS"),
    desc_block(data, father_pgs_vars, rep(NA_character_, length(father_pgs_vars)), father_pgs_labels, "Father PGS"))
}

trio_sample <- df %>% filter(!is.na(PC1) &  !is.na(PC1_fp) &  !is.na(PC1_mp))
duo_sample  <- df %>% filter(!is.na(PC1) & (!is.na(PC1_fp) | !is.na(PC1_mp)))

trio_desc <- desc_sample(trio_sample)
duo_desc  <- desc_sample(duo_sample)

pgs_family_cor_trio_num <- family_cor(trio_sample)
write.csv(pgs_family_cor_trio_num, file.path(outFileStem, "GN_PGS_family_correlations_trios.csv"), row.names = FALSE)

pgs_family_cor_trio <- fmt_family_cor(family_cor(trio_sample))
wb_fam <- createWorkbook()
addWorksheet(wb_fam, "PGS family corr (trios)")
writeData(wb_fam, "PGS family corr (trios)", pgs_family_cor_trio)
saveWorkbook(wb_fam, file.path(outFileStem, "GN_PGS_family_correlations_trios_maintext.xlsx"),
             overwrite = TRUE)

g1 <- "Child + >=1 Parent DNA"; g2 <- "Child Only"; g3 <- "No Child DNA"

df_g1g2 <- df %>%
  mutate(DNA_Group = case_when(
    !is.na(PC1) & (!is.na(PC1_fp) | !is.na(PC1_mp)) ~ g1,
    !is.na(PC1) &   is.na(PC1_fp) &  is.na(PC1_mp)  ~ g2,
    TRUE ~ NA_character_)) %>% filter(!is.na(DNA_Group))
fam_with_any_dna <- df_full %>%
  filter(!is.na(PC1) | !is.na(PC1_fp) | !is.na(PC1_mp)) %>% distinct(id_fam)
df_g3 <- df_full %>% anti_join(fam_with_any_dna, by = "id_fam") %>%
  group_by(id_fam) %>% arrange(random, .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
  mutate(DNA_Group = g3)
stopifnot(length(intersect(df_g3$id_fam, df_g1g2$id_fam)) == 0)
df_all <- bind_rows(df_g1g2, df_g3) %>%
  mutate(DNA_Group = factor(DNA_Group, levels = c(g1, g2, g3)))

compare_pair <- function(data, vars, ga, gb) {
  map_dfr(intersect(vars, names(data)), function(v) {
    a <- data %>% filter(DNA_Group == ga) %>% pull(!!sym(v)) %>% na.omit()
    b <- data %>% filter(DNA_Group == gb) %>% pull(!!sym(v)) %>% na.omit()
    if (length(a) < 2 || length(b) < 2) return(NULL)
    tt  <- t.test(a, b)
    psd <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
    d   <- if (isTRUE(psd > 0)) (mean(a) - mean(b)) / psd else NA_real_
    tibble(Variable = v, Group_A = ga, Group_B = gb,
           N_A = length(a), Mean_A = mean(a), SD_A = sd(a),
           Skew_A = if (length(a) > 2) skewness(a) else NA_real_,
           N_B = length(b), Mean_B = mean(b), SD_B = sd(b),
           Skew_B = if (length(b) > 2) skewness(b) else NA_real_,
           Cohen_d = d,
           t_stat = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value)
  }) %>%
    mutate(p_adjusted_FDR = p.adjust(p_value, "fdr"),
           Trait_Type = var_type_map[Variable], Label = var_label_map[Variable]) %>%
    relocate(Trait_Type, Label, .before = Variable)
}

cmp_G1_G2 <- compare_pair(df_all, c(pheno_vars, ses_vars, child_pgs_vars), g1, g2)
cmp_G2_G3 <- compare_pair(df_all, c(pheno_vars, ses_vars),                 g2, g3)

write.csv(trio_desc, file.path(outFileStem, "GN_Unrelated_Trios_descriptives.csv"), row.names = FALSE)
write.csv(duo_desc,  file.path(outFileStem, "GN_Unrelated_allDuos_descriptives.csv"),  row.names = FALSE)
write.csv(cmp_G1_G2, file.path(outFileStem, "GN_Compare_G1wParentDNA_G2noParentDNA.csv"), row.names = FALSE)
write.csv(cmp_G2_G3, file.path(outFileStem, "GN_Compare_G2noParentDNA_G3noDNA.csv"), row.names = FALSE)

sample_char <- function(data, label) {
  tibble(
    `Sample (unrelated: one child per family)` = label,
    N           = nrow(data),
    Prop_Female = if ("sex1"  %in% names(data)) mean(data$sex1  == 0, na.rm = TRUE) else NA_real_,
    Prop_DZ     = if ("zygos" %in% names(data)) mean(data$zygos == 2, na.rm = TRUE) else NA_real_)
}
sample_characteristics <- bind_rows(
  sample_char(df_all,                             "Whole sample"),
  sample_char(trio_sample,                        "Completed trios (child+mother+father)"),
  sample_char(duo_sample,                         "Duos (child + >=1 parent)"),
  sample_char(df_all %>% filter(DNA_Group == g2), "Child DNA only"),
  sample_char(df_all %>% filter(DNA_Group == g3), "No DNA"))
write.csv(sample_characteristics, file.path(outFileStem, "GN_Sample_characteristics.csv"), row.names = FALSE)

sensi_group_covars <- intersect(c("sex1", "zygos", "chiptype"), names(df))
sensi_cont_covars  <- character(0)

cov_disp <- c(sex1 = "sex", zygos = "zygosity", chiptype = "genotyping chip", age = "age")

sensi_one <- function(y, x, data, x_factor) {
  d <- data[!is.na(data[[y]]) & !is.na(data[[x]]), ]
  if (nrow(d) < 10) return(NULL)
  rhs <- if (x_factor) sprintf("factor(%s)", x) else x
  m <- lm(as.formula(paste(y, "~", rhs)), data = d)
  f <- summary(m)$fstatistic
  if (is.null(f)) return(NULL)
  tibble(Variable = y, N = nrow(d),
         F_stat = unname(f[1]),
         P_Value = pf(f[1], f[2], f[3], lower.tail = FALSE),
         R2 = summary(m)$r.squared)
}

sensi_targets <- bind_rows(
  tibble(var = pheno_vars, type = pheno_types, label = pheno_labels,
         age_col = sapply(VarPairs, function(x) x[3])),
  tibble(var = ses_vars, type = "SES", label = ses_labels, age_col = NA_character_))

sensi_rows <- list()
for (ti in seq_len(nrow(sensi_targets))) {
  y <- sensi_targets$var[ti]
  if (!y %in% names(df)) next
  age_col <- sensi_targets$age_col[ti]
  specs <- c(setNames(rep(TRUE,  length(sensi_group_covars)), sensi_group_covars),
             setNames(rep(FALSE, length(sensi_cont_covars)),  sensi_cont_covars))
  if (!is.na(age_col) && age_col %in% names(df)) specs <- c(specs, setNames(FALSE, age_col))
  for (cv in names(specs)) {
    r <- sensi_one(y, cv, df, specs[[cv]])
    if (is.null(r)) next
    r$Covariate  <- if (!is.na(age_col) && cv == age_col) "age" else cv
    r$Trait_Type <- sensi_targets$type[ti]; r$Label <- sensi_targets$label[ti]
    sensi_rows[[length(sensi_rows) + 1]] <- r
  }
}
sensi <- bind_rows(sensi_rows) %>%
  group_by(Covariate) %>% mutate(P_adjusted_FDR = p.adjust(P_Value, "fdr")) %>% ungroup() %>%
  mutate(Covariate = dplyr::recode(Covariate, !!!cov_disp)) %>%
  select(Covariate, Trait_Type, Label, Variable, N, F_stat, P_Value, P_adjusted_FDR, R2)

write.csv(sensi, file.path(outFileStem, "GN_Sensitivity.csv"), row.names = FALSE)

pgs_all    <- c(child_pgs_vars, mother_pgs_vars, father_pgs_vars)
pgs_member <- rep(c("Child", "Mother", "Father"), each = length(child_pgs_vars))
pgs_trait  <- rep(unname(pgs_short[child_pgs_vars]), times = 3)

pgs_sensi_rows <- list()
for (k in seq_along(pgs_all)) {
  y <- pgs_all[k]
  if (!y %in% names(df)) next
  specs <- c(setNames(rep(TRUE,  length(sensi_group_covars)), sensi_group_covars),
             setNames(rep(FALSE, length(sensi_cont_covars)),  sensi_cont_covars))
  for (cv in names(specs)) {
    r <- sensi_one(y, cv, df, specs[[cv]])
    if (is.null(r)) next
    r$Covariate <- cv; r$Member <- pgs_member[k]; r$Trait <- pgs_trait[k]
    pgs_sensi_rows[[length(pgs_sensi_rows) + 1]] <- r
  }
}
pgs_sensi <- bind_rows(pgs_sensi_rows) %>%
  group_by(Covariate) %>% mutate(P_adjusted_FDR = p.adjust(P_Value, "fdr")) %>% ungroup() %>%
  mutate(Covariate = dplyr::recode(Covariate, !!!cov_disp)) %>%
  select(Covariate, Member, Trait, Variable, N, F_stat, P_Value, P_adjusted_FDR, R2)

write.csv(pgs_sensi, file.path(outFileStem, "GN_PGS_Sensitivity.csv"), row.names = FALSE)

var_to_label <- setNames(corr_map$Label, corr_map$Var)
ses_lab      <- var_to_label[ses_vars]; ses_lab <- ses_lab[!is.na(ses_lab)]
ses_time_lab <- c("SES (Birth)" = "Birth", "SES (Age 7)" = "Age 7",
                  "SES (Age 16)" = "Age 16", "SES (Age 21)" = "Age 21")

pgs_rows_var <- c(child_pgs_vars, mother_pgs_vars, father_pgs_vars)
pgs_rows_var <- pgs_rows_var[pgs_rows_var %in% names(var_to_label)]
pgs_rows_lab <- var_to_label[pgs_rows_var]

pgs_ses_df <- as.data.frame(corr_matrix[pgs_rows_lab, ses_lab, drop = FALSE]) %>%
  rownames_to_column("Label") %>%
  pivot_longer(-Label, names_to = "SES_Timepoint", values_to = "r") %>%
  mutate(Member = case_when(grepl("^PGS child ",  Label) ~ "Child",
                            grepl("^PGS mother ", Label) ~ "Mother",
                            grepl("^PGS father ", Label) ~ "Father"),
         Member = factor(Member, levels = c("Child", "Mother", "Father")),
         Trait  = sub("^PGS (child|mother|father) \\([^)]*\\): ", "", Label),
         Trait  = factor(Trait, levels = unname(pgs_short[child_pgs_vars])),
         SES_Timepoint = factor(unname(ses_time_lab[SES_Timepoint]), levels = unname(ses_time_lab)))

p_pgs_ses <- ggplot(pgs_ses_df, aes(SES_Timepoint, Trait, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))), size = 4) +
  facet_wrap(~ Member, ncol = 3) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey90", name = "r") +
  labs(title = "PGS-SES correlations by family member", x = NULL, y = NULL) +
  theme_gn_corr()
save_gn(p_pgs_ses, "SES_Corr_PGS", width = 14, height = 8)

pheno_domain_map <- setNames(sapply(VarPairs, function(x) x[5]), pheno_vars)

pheno_rows_var <- pheno_vars[pheno_vars %in% names(var_to_label)]
pheno_rows_lab <- var_to_label[pheno_rows_var]

pheno_ses_df <- as.data.frame(corr_matrix[pheno_rows_lab, ses_lab, drop = FALSE]) %>%
  rownames_to_column("Label") %>%
  pivot_longer(-Label, names_to = "SES_Timepoint", values_to = "r") %>%
  left_join(tibble(Label = pheno_rows_lab, Variable = pheno_rows_var,
                   Domain = unname(pheno_domain_map[pheno_rows_var]),
                   PlainLabel = unname(var_label_map[pheno_rows_var])),
            by = "Label") %>%
  mutate(SES_Timepoint = factor(unname(ses_time_lab[SES_Timepoint]), levels = unname(ses_time_lab)),
         PlainLabel = factor(PlainLabel, levels = unique(unname(var_label_map[pheno_rows_var]))))

p_pheno_ses <- ggplot(pheno_ses_df, aes(SES_Timepoint, PlainLabel, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))), size = 6) +
  facet_wrap(~ Domain, scales = "free_y") +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey90", name = "r") +
  labs(title = "Child phenotype-SES correlations by trait domain", x = NULL, y = NULL) +
  theme_gn_corr()
save_gn(p_pheno_ses, "SES_Corr_Phenotypes", width = 16, height = 20)

wb_all <- createWorkbook()
add_sheet <- function(name, x, rn = FALSE) {
  addWorksheet(wb_all, name); writeData(wb_all, name, x, rowNames = rn)
}
add_sheet("Correlation matrix", corr_matrix, rn = TRUE)
freezePane(wb_all, "Correlation matrix", firstActiveRow = 2, firstActiveCol = 2)
add_sheet("Unrelated Trios (Complete)", trio_desc)
add_sheet("Unrelated Duos (>=1 parent)", duo_desc)
add_sheet("G1wParentDNA vs G2noParentDNA", cmp_G1_G2)
add_sheet("G2noParentDNA vs G3noDNA",      cmp_G2_G3)
add_sheet("Sample characteristics", sample_characteristics)
add_sheet("Pheno+SES sensitivity", sensi)
add_sheet("PGS sensitivity", pgs_sensi)
saveWorkbook(wb_all, file.path(outFileStem, "GN_DataChecking.xlsx"), overwrite = TRUE)
