library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)

source("0GN_Theme.R")

outFileStem <- "GN_Results/Reg/"

GN_Reg <- read.csv(paste0(outFileStem, "GN_Reg.csv")) %>%
  filter(Subgroup == "All_Trios")

trait_order <- unique(GN_Reg$Unique_Label)

domain_order <- c(
  "Cognitive Abilities",
  "Education Achievement & Attainment",
  "Height & BMI",
  "Psychopathology"
)

lighten <- function(hex, amt = 0.55) lighten_gn(hex, f = amt)

fam_cor <- read.csv("GN_Results/Data_Checking/GN_PGS_family_correlations_trios.csv") %>%
  filter(Trait != "N (pairs)") %>%
  transmute(pgs_label = Trait, r_mc = r_MotherChild, r_fc = r_FatherChild)

source("0GN_VarPairs.R")
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
trait_pgs <- data.frame(trait     = sapply(VarPairs, function(x) x[6]),
                        pgs_label = unname(pgs_short[sapply(VarPairs, function(x) x[2])]))

dec <- GN_Reg %>%
  transmute(trait  = Unique_Label,
            Domain = factor(Domain, levels = domain_order),
            pop    = Beta_Child_Pop_Est,
            direct = Beta_Child_OverallGN_Est,
            mom    = Beta_Mom_OverallGN_r_Est,
            dad    = Beta_Dad_OverallGN_r_Est) %>%
  left_join(trait_pgs, by = "trait") %>%
  left_join(fam_cor, by = "pgs_label") %>%
  mutate(mom_w = r_mc * mom,
         dad_w = r_fc * dad)

LAB_DIR   <- "Child PGS (Direct)"
LAB_MOM_W <- "Maternal PGS (Genetic Nurture), weighted by mother-child PGS correlation"
LAB_DAD_W <- "Paternal PGS (Genetic Nurture), weighted by father-child PGS correlation"
LAB_POP   <- "Child PGS (Population)"
BAR_COL_W <- setNames(c(CHD_DIR_COL, MAT_COL, PAT_COL), c(LAB_DIR, LAB_MOM_W, LAB_DAD_W))

dec <- dec %>%
  mutate(trait_clean = paste0(toupper(substr(trait, 1, 1)), substr(trait, 2, nchar(trait)))) %>%
  group_by(Domain) %>%
  mutate(idx = match(trait, rev(trait_order[trait_order %in% trait]))) %>%
  ungroup()
dec_lvl <- unique(dec$trait_clean[order(match(dec$trait, rev(trait_order)))])
dec$trait_clean <- factor(dec$trait_clean, levels = dec_lvl)

H        <- 0.5
POP_LN_H <- 0.7
dec_seg <- dec %>%
  select(trait, trait_clean, Domain, idx, direct, mom_w, dad_w) %>%
  pivot_longer(c(direct, mom_w, dad_w), names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(comp, levels = c("direct", "mom_w", "dad_w"))) %>%
  arrange(trait, comp) %>%
  group_by(trait) %>%
  mutate(pos = pmax(val, 0), neg = pmin(val, 0),
         xmax = ifelse(val >= 0, cumsum(pos), cumsum(neg) - neg),
         xmin = ifelse(val >= 0, cumsum(pos) - pos, cumsum(neg))) %>%
  ungroup() %>%
  mutate(series = factor(recode(comp, direct = LAB_DIR, mom_w = LAB_MOM_W, dad_w = LAB_DAD_W),
                         levels = names(BAR_COL_W)),
         ymin = idx - H / 2, ymax = idx + H / 2)
dec_popline <- dec %>%
  transmute(Domain, trait_clean, x = pop,
            ymin = idx - POP_LN_H / 2, ymax = idx + POP_LN_H / 2)

plot_main <- ggplot(dec_seg) +
  geom_blank(aes(y = trait_clean)) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey50") +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = series),
            colour = "white", linewidth = 0.2) +
  geom_segment(data = dec_popline, aes(x = x, xend = x, y = ymin, yend = ymax, linetype = LAB_POP),
               colour = "black", linewidth = 0.8, key_glyph = "vpath") +
  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = BAR_COL_W, name = NULL) +
  scale_linetype_manual(values = setNames("solid", LAB_POP), name = NULL) +
  scale_x_continuous(breaks = BRK_BETA) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  labs(x = expression("Standardised" ~ beta), y = NULL,
       title = "Population prediction and its decomposition into direct and genetic nurture effects") +
  theme_gn(grid = "x") +
  guides(fill = guide_legend(ncol = 1, order = 1), linetype = guide_legend(order = 2))

dec_height <- 4.0 + ((length(trait_order) / 2) * 0.8)
save_gn(plot_main, "Reg_PopDecomposition_Bar", width = 24, height = dec_height)

construct_map <- c(
  "primary school grades age 7"  = "Educational achievement",
  "primary school grades age 9"  = "Educational achievement",
  "primary school grades age 10" = "Educational achievement",
  "primary school grades age 12" = "Educational achievement",
  "GCSE grades age 16"           = "Educational achievement",
  "A-level grades age 18"        = "Educational achievement",
  "university grades age 21"     = "Educational achievement",
  "A-level enrolment age 16"    = "Educational attainment",
  "university enrolment age 18" = "Educational attainment",
  "years of schooling age 26"    = "Educational attainment",
  "g age 7"  = "General cognitive ability (g)", "g age 9"  = "General cognitive ability (g)",
  "g age 10" = "General cognitive ability (g)", "g age 12" = "General cognitive ability (g)",
  "g age 16" = "General cognitive ability (g)", "g age 25" = "General cognitive ability (g)",
  "verbal ability age 7"  = "Verbal ability", "verbal ability age 9"  = "Verbal ability",
  "verbal ability age 10" = "Verbal ability", "verbal ability age 12" = "Verbal ability",
  "verbal ability age 16" = "Verbal ability", "verbal ability age 25" = "Verbal ability",
  "nonverbal ability age 7"  = "Nonverbal ability", "nonverbal ability age 9"  = "Nonverbal ability",
  "nonverbal ability age 10" = "Nonverbal ability", "nonverbal ability age 12" = "Nonverbal ability",
  "nonverbal ability age 16" = "Nonverbal ability", "nonverbal ability age 25" = "Nonverbal ability",
  "height age 7"  = "Height", "height age 12" = "Height", "height age 14" = "Height",
  "height age 16" = "Height", "height age 22" = "Height", "height age 26" = "Height",
  "BMI age 7"  = "BMI", "BMI age 12" = "BMI", "BMI age 14" = "BMI",
  "BMI age 16" = "BMI", "BMI age 22" = "BMI", "BMI age 26" = "BMI"
)
psy_traits <- dec$trait[dec$Domain == "Psychopathology"]
construct_map <- c(construct_map, setNames(psy_traits, psy_traits))

dec_w <- dec %>%
  select(trait, Domain, pop, direct, mom_w, dad_w)

construct_order_w <- c("Educational attainment", "Educational achievement",
                       "General cognitive ability (g)", "Verbal ability", "Nonverbal ability",
                       "Psychopathology", "Height", "BMI")

dec_w_agg <- dec_w %>%
  mutate(construct = ifelse(as.character(Domain) == "Psychopathology", "Psychopathology",
                            unname(construct_map[trait]))) %>%
  group_by(construct) %>%
  summarise(n_traits = n(),
            pop    = mean(pop),    direct = mean(direct),
            mom_w  = mean(mom_w),  dad_w  = mean(dad_w),
            .groups = "drop")

LAB_MOM_W <- "Maternal PGS (Genetic Nurture), weighted by mother-child PGS correlation"
LAB_DAD_W <- "Paternal PGS (Genetic Nurture), weighted by father-child PGS correlation"
BAR_COL_W <- setNames(c(CHD_DIR_COL, MAT_COL, PAT_COL), c(LAB_DIR, LAB_MOM_W, LAB_DAD_W))
H_W      <- 0.5
POP_LN_H <- 0.7

dec_w_agg <- dec_w_agg %>%
  mutate(trait = construct,
         trait_clean = paste0(toupper(substr(trait, 1, 1)), substr(trait, 2, nchar(trait))),
         idx = match(trait, rev(construct_order_w)))
w_lvl <- unique(dec_w_agg$trait_clean[order(match(dec_w_agg$trait, rev(construct_order_w)))])
dec_w_agg$trait_clean <- factor(dec_w_agg$trait_clean, levels = w_lvl)

w_seg <- dec_w_agg %>%
  select(trait, trait_clean, idx, direct, mom_w, dad_w) %>%
  pivot_longer(c(direct, mom_w, dad_w), names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(comp, levels = c("direct", "mom_w", "dad_w"))) %>%
  arrange(trait, comp) %>%
  group_by(trait) %>%
  mutate(pos = pmax(val, 0), neg = pmin(val, 0),
         xmax = ifelse(val >= 0, cumsum(pos), cumsum(neg) - neg),
         xmin = ifelse(val >= 0, cumsum(pos) - pos, cumsum(neg))) %>%
  ungroup() %>%
  mutate(series = factor(recode(comp, direct = LAB_DIR, mom_w = LAB_MOM_W, dad_w = LAB_DAD_W),
                         levels = names(BAR_COL_W)),
         ymin = idx - H_W / 2, ymax = idx + H_W / 2)
w_pop <- dec_w_agg %>%
  transmute(trait_clean, x = pop,
            ymin = idx - POP_LN_H / 2, ymax = idx + POP_LN_H / 2)

plot_main_agg_w <- ggplot(w_seg) +
  geom_blank(aes(y = trait_clean)) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey50") +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = series),
            colour = "white", linewidth = 0.2) +
  geom_segment(data = w_pop, aes(x = x, xend = x, y = ymin, yend = ymax, linetype = LAB_POP),
               colour = "black", linewidth = 0.8, key_glyph = "vpath") +
  scale_fill_manual(values = BAR_COL_W, name = NULL) +
  scale_linetype_manual(values = setNames("solid", LAB_POP), name = NULL) +
  scale_x_continuous(breaks = BRK_BETA) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  labs(x = expression("Standardised" ~ beta), y = NULL,
       title = "Population prediction and its decomposition into direct and genetic nurture effects, averaged across ages") +
  theme_gn(grid = "x") +
  guides(fill = guide_legend(ncol = 1, order = 1), linetype = guide_legend(order = 2))

w_height <- 3.0 + nrow(dec_w_agg) * 0.55
save_gn(plot_main_agg_w, "Reg_PopDecomposition_Bar_Agg_rWeighted", width = 16, height = w_height)

forest_colors <- SERIES4_COL

p2_data <- GN_Reg %>%
  select(Unique_Label, Domain,
         est_pop   = Beta_Child_Pop_Est,
         est_child = Beta_Child_OverallGN_Est,
         est_mat   = Beta_Mom_OverallGN_r_Est,
         est_pat   = Beta_Dad_OverallGN_r_Est,

         lo_pop    = Beta_Child_Pop_LowerCI,
         hi_pop    = Beta_Child_Pop_UpperCI,
         lo_child  = Beta_Child_OverallGN_LowerCI,
         hi_child  = Beta_Child_OverallGN_UpperCI,
         lo_mat    = Beta_Mom_OverallGN_r_LowerCI,
         hi_mat    = Beta_Mom_OverallGN_r_UpperCI,
         lo_pat    = Beta_Dad_OverallGN_r_LowerCI,
         hi_pat    = Beta_Dad_OverallGN_r_UpperCI,

         p_pop     = Beta_Child_Pop_P_FDR,
         p_child   = Beta_Child_OverallGN_P_FDR,
         p_mat     = Beta_Mom_OverallGN_P_FDR,
         p_pat     = Beta_Dad_OverallGN_P_FDR
  ) %>%
  pivot_longer(
    cols      = -c(Unique_Label, Domain),
    names_to  = c(".value", "Model"),
    names_pattern = "(est|lo|hi|p)_(.*)"
  ) %>%
  mutate(
    Model = dplyr::recode(Model,
                          pop   = "Child PGS (Population)",
                          child = "Child PGS (Direct)",
                          mat   = "Maternal PGS (Genetic Nurture)",
                          pat   = "Paternal PGS (Genetic Nurture)"),

    Model = factor(Model, levels = c("Paternal PGS (Genetic Nurture)", "Maternal PGS (Genetic Nurture)", "Child PGS (Direct)", "Child PGS (Population)")),

    trait_clean = paste0(toupper(substr(Unique_Label, 1, 1)), substr(Unique_Label, 2, nchar(Unique_Label))),

    Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
    trait_clean  = factor(trait_clean,
                          levels = unique(trait_clean[order(match(Unique_Label, levels(Unique_Label)))])),
    Domain       = factor(Domain, levels = domain_order),

    Stars = fdr_stars(p)
  )

dodge_width <- 0.65

plot2 <- ggplot(p2_data, aes(x = est, y = trait_clean, color = Model)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey50") +

  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = dodge_width), linewidth = 0.8) +

  geom_point(size = 3, position = position_dodge(width = dodge_width)) +

  geom_text(aes(x = hi + max(abs(hi), na.rm = TRUE) * 0.05, label = Stars),
            position = position_dodge(width = dodge_width),
            vjust = 0.75, size = 4.5, fontface = "bold", show.legend = FALSE) +

  scale_color_manual(values = forest_colors,
                     breaks = c("Child PGS (Population)", "Child PGS (Direct)", "Maternal PGS (Genetic Nurture)", "Paternal PGS (Genetic Nurture)")) +
  scale_x_continuous(breaks = BRK_BETA) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +

  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +

  labs(
    title    = "Child PGS effects and parental genetic nurture effects",
    x        = X_BETA,
    y        = NULL,
    color    = NULL
  ) +
  theme_gn(grid = "x") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

calc_height <- 4.0 + ((length(unique(p2_data$Unique_Label)) / 2) * 0.8)

save_gn(plot2, "Reg_PGSDecomposition", width = 24, height = calc_height)

gn_colors <- c(
  "Maternal Genetic Nurture (Duo vs Population Baseline)" = lighten(MAT_COL),
  "Maternal Genetic Nurture (Trio vs Duo)"     = MAT_COL,
  "Paternal Genetic Nurture (Duo vs Population Baseline)" = lighten(PAT_COL),
  "Paternal Genetic Nurture (Trio vs Duo)"     = PAT_COL,
  "Overall Genetic Nurture"                   = OVR_COL
)

p3_data <- GN_Reg %>%
  select(Unique_Label, Domain,
         `Maternal Genetic Nurture (Duo vs Population Baseline)` = DeltaR2_Mat_Marginal_Est,
         `Maternal Genetic Nurture (Trio vs Duo)`     = DeltaR2_Unique_Mat_Est,
         `Paternal Genetic Nurture (Duo vs Population Baseline)` = DeltaR2_Pat_Marginal_Est,
         `Paternal Genetic Nurture (Trio vs Duo)`     = DeltaR2_Unique_Pat_Est,
         `Overall Genetic Nurture`        = DeltaR2_Overall_GN_Est,
         p_mat_marg = WaldF_Mat_Marginal_P,
         p_mat_uniq = WaldF_Unique_Mat_P,
         p_pat_marg = WaldF_Pat_Marginal_P,
         p_pat_uniq = WaldF_Unique_Pat_P,
         p_overall  = WaldF_Overall_P) %>%
  pivot_longer(cols = c(`Maternal Genetic Nurture (Duo vs Population Baseline)`, `Maternal Genetic Nurture (Trio vs Duo)`,
                        `Paternal Genetic Nurture (Duo vs Population Baseline)`, `Paternal Genetic Nurture (Trio vs Duo)`, `Overall Genetic Nurture`),
               names_to = "Series", values_to = "DeltaR2") %>%
  mutate(
    P_Value = case_when(
      Series == "Maternal Genetic Nurture (Duo vs Population Baseline)" ~ p_mat_marg,
      Series == "Maternal Genetic Nurture (Trio vs Duo)"      ~ p_mat_uniq,
      Series == "Paternal Genetic Nurture (Duo vs Population Baseline)"  ~ p_pat_marg,
      Series == "Paternal Genetic Nurture (Trio vs Duo)"      ~ p_pat_uniq,
      Series == "Overall Genetic Nurture"         ~ p_overall,
      TRUE                            ~ NA_real_
    ),
    Stars = fdr_stars(P_Value),
    Series = factor(Series, levels = c("Overall Genetic Nurture",
                                       "Paternal Genetic Nurture (Duo vs Population Baseline)", "Paternal Genetic Nurture (Trio vs Duo)",
                                       "Maternal Genetic Nurture (Duo vs Population Baseline)", "Maternal Genetic Nurture (Trio vs Duo)")),
    Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
    Domain       = factor(Domain, levels = domain_order)
  )

plot3 <- ggplot(p3_data, aes(x = DeltaR2, y = Unique_Label, fill = Series)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_col(position = position_dodge(width = 0.75), width = 0.72) +
  geom_text(aes(x = DeltaR2 + max(abs(DeltaR2), na.rm = TRUE) * 0.04, label = Stars),
            position = position_dodge(width = 0.75),
            size = 4.5, hjust = 0, fontface = "bold") +
  scale_fill_manual(
    values = gn_colors,
    guide  = guide_legend(reverse = TRUE)
  ) +
  scale_y_discrete(
    labels = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  ) +
  scale_x_continuous(breaks = BRK_R2, expand = expansion(mult = c(0.05, 0.15)),
                     labels = scales::percent_format(accuracy = 0.5)) +
  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +
  labs(
    title = "Omitted-parent effects: incremental R² for genetic nurture effects from duo versus trio models",
    x    = expression(Delta*R^2),
    y    = NULL,
    fill = ""
  ) +
  theme_gn(grid = "x") +
  guides(fill = guide_legend(reverse = TRUE, nrow = 3, byrow = TRUE)) +
  coord_cartesian(clip = "off")

calc_height_p3 <- 24

save_gn(plot3, "RegSup1_deltaR2", width = 24, height = calc_height_p3)

base_cols <- c("Child PGS"     = CHD_DIR_COL,
               "Maternal PGS"  = MAT_COL,
               "Paternal PGS"  = PAT_COL)

p4_levels <- c("Paternal PGS (noPC Trios)",  "Paternal PGS (All Trios)",
               "Maternal PGS (noPC Trios)",  "Maternal PGS (All Trios)",
               "Child PGS (noPC Trios)",     "Child PGS (All Trios)")
p4_colors <- setNames(
  c(lighten(base_cols[["Paternal PGS"]]), base_cols[["Paternal PGS"]],
    lighten(base_cols[["Maternal PGS"]]), base_cols[["Maternal PGS"]],
    lighten(base_cols[["Child PGS"]]),    base_cols[["Child PGS"]]),
  p4_levels)

GN_Reg_pccmp <- read.csv(paste0(outFileStem, "GN_Reg.csv")) %>%
  filter(Subgroup %in% c("All_Trios", "noPC_Trios"))

p4_data <- GN_Reg_pccmp %>%
  select(Unique_Label, Domain, Subgroup,
         est_child = Beta_Child_OverallGN_Est,
         est_mat   = Beta_Mom_OverallGN_r_Est,
         est_pat   = Beta_Dad_OverallGN_r_Est,
         lo_child  = Beta_Child_OverallGN_LowerCI,
         hi_child  = Beta_Child_OverallGN_UpperCI,
         lo_mat    = Beta_Mom_OverallGN_r_LowerCI,
         hi_mat    = Beta_Mom_OverallGN_r_UpperCI,
         lo_pat    = Beta_Dad_OverallGN_r_LowerCI,
         hi_pat    = Beta_Dad_OverallGN_r_UpperCI,
         p_child   = Beta_Child_OverallGN_P_FDR,
         p_mat     = Beta_Mom_OverallGN_P_FDR,
         p_pat     = Beta_Dad_OverallGN_P_FDR) %>%
  pivot_longer(cols = -c(Unique_Label, Domain, Subgroup),
               names_to = c(".value", "EstKey"),
               names_pattern = "(est|lo|hi|p)_(.*)") %>%
  mutate(
    Estimate = dplyr::recode(EstKey,
                             child = "Child PGS",
                             mat   = "Maternal PGS",
                             pat   = "Paternal PGS"),
    Sample = ifelse(Subgroup == "All_Trios", "All Trios", "noPC Trios"),
    Group  = paste0(Estimate, " (", Sample, ")"),
    Group = factor(Group, levels = p4_levels),
    Stars = fdr_stars(p),
    Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
    Domain       = factor(Domain, levels = domain_order)
  )

dodge_width_p4 <- 0.7

plot4 <- ggplot(p4_data, aes(x = est, y = Unique_Label, color = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = dodge_width_p4), linewidth = 0.8) +
  geom_point(size = 2.6, position = position_dodge(width = dodge_width_p4)) +
  geom_text(aes(x = hi + max(abs(hi), na.rm = TRUE) * 0.05, label = Stars),
            position = position_dodge(width = dodge_width_p4),
            vjust = 0.75, size = 4.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = p4_colors, guide = guide_legend(reverse = TRUE, nrow = 3, byrow = TRUE)) +
  scale_y_discrete(
    labels = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  ) +
  scale_x_continuous(breaks = BRK_BETA, expand = expansion(mult = c(0.05, 0.10))) +
  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +
  labs(x = X_BETA, y = NULL, color = "",
       title = "Sensitivity to population structure: direct child PGS and genetic nurture effects, with versus without principal components") +
  theme_gn(grid = "x") +
  coord_cartesian(clip = "off")

calc_height_p4 <- 28

save_gn(plot4, "RegSup2_PC", width = 24, height = calc_height_p4)

p5_levels <- c("Paternal PGS (Duo)",  "Paternal PGS (Trio)",
               "Maternal PGS (Duo)",  "Maternal PGS (Trio)")
p5_colors <- setNames(
  c(lighten(PAT_COL), PAT_COL, lighten(MAT_COL), MAT_COL),
  p5_levels
)

p5_data <- GN_Reg %>%
  select(Unique_Label, Domain,
         est_matduo  = Beta_Mom_MatGN_r_Est,
         est_mattrio = Beta_Mom_OverallGN_r_Est,
         est_patduo  = Beta_Dad_PatGN_r_Est,
         est_pattrio = Beta_Dad_OverallGN_r_Est,
         lo_matduo   = Beta_Mom_MatGN_r_LowerCI,
         hi_matduo   = Beta_Mom_MatGN_r_UpperCI,
         lo_mattrio  = Beta_Mom_OverallGN_r_LowerCI,
         hi_mattrio  = Beta_Mom_OverallGN_r_UpperCI,
         lo_patduo   = Beta_Dad_PatGN_r_LowerCI,
         hi_patduo   = Beta_Dad_PatGN_r_UpperCI,
         lo_pattrio  = Beta_Dad_OverallGN_r_LowerCI,
         hi_pattrio  = Beta_Dad_OverallGN_r_UpperCI,
         p_matduo    = Beta_Mom_MatGN_P_FDR,
         p_mattrio   = Beta_Mom_OverallGN_P_FDR,
         p_patduo    = Beta_Dad_PatGN_P_FDR,
         p_pattrio   = Beta_Dad_OverallGN_P_FDR) %>%
  pivot_longer(cols = -c(Unique_Label, Domain),
               names_to = c(".value", "Key"),
               names_pattern = "(est|lo|hi|p)_(.*)") %>%
  mutate(
    Series = dplyr::recode(Key,
                           matduo  = "Maternal PGS (Duo)",
                           mattrio = "Maternal PGS (Trio)",
                           patduo  = "Paternal PGS (Duo)",
                           pattrio = "Paternal PGS (Trio)"),
    Series = factor(Series, levels = p5_levels),
    Stars = fdr_stars(p),
    Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
    Domain       = factor(Domain, levels = domain_order)
  )

dodge_width_p5 <- 0.7

plot5 <- ggplot(p5_data, aes(x = est, y = Unique_Label, color = Series)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = dodge_width_p5), linewidth = 0.8) +
  geom_point(size = 2.8, position = position_dodge(width = dodge_width_p5)) +
  geom_text(aes(x = hi + max(abs(hi), na.rm = TRUE) * 0.05, label = Stars),
            position = position_dodge(width = dodge_width_p5),
            vjust = 0.75, size = 4.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = p5_colors,
                     guide = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE)) +
  scale_y_discrete(
    labels = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  ) +
  scale_x_continuous(breaks = BRK_BETA, expand = expansion(mult = c(0.05, 0.10))) +
  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +
  labs(x = X_BETA, y = NULL, color = "",
       title = "Omitted-parent effects: beta estimates for genetic nurture effects from duo versus trio models") +
  theme_gn(grid = "x") +
  coord_cartesian(clip = "off")

save_gn(plot5, "RegSup3_TrioDuoBeta", width = 24, height = 24)

p6_levels <- c("Paternal PGS (Duos in Complete Trios)", "Paternal PGS (All Duos)",
               "Maternal PGS (Duos in Complete Trios)", "Maternal PGS (All Duos)")
p6_colors <- setNames(
  c(lighten(PAT_COL), PAT_COL, lighten(MAT_COL), MAT_COL),
  p6_levels
)

GN_Reg_duofull <- read.csv(paste0(outFileStem, "GN_Reg.csv")) %>%
  filter(Subgroup %in% c("All_Trios", "Maternal_Duo_Full", "Paternal_Duo_Full"))

pull_series <- function(subgroup, est_col, lo_col, hi_col, p_col, series_name) {
  GN_Reg_duofull %>%
    filter(Subgroup == subgroup) %>%
    transmute(Unique_Label, Domain,
              est = .data[[est_col]], lo = .data[[lo_col]], hi = .data[[hi_col]], p = .data[[p_col]],
              Series = series_name)
}

p6_data <- bind_rows(
  pull_series("All_Trios",        "Beta_Mom_MatGN_r_Est", "Beta_Mom_MatGN_r_LowerCI", "Beta_Mom_MatGN_r_UpperCI", "Beta_Mom_MatGN_P_FDR", "Maternal PGS (Duos in Complete Trios)"),
  pull_series("Maternal_Duo_Full",   "Beta_Mom_MatGN_r_Est", "Beta_Mom_MatGN_r_LowerCI", "Beta_Mom_MatGN_r_UpperCI", "Beta_Mom_MatGN_P_FDR", "Maternal PGS (All Duos)"),
  pull_series("All_Trios",        "Beta_Dad_PatGN_r_Est", "Beta_Dad_PatGN_r_LowerCI", "Beta_Dad_PatGN_r_UpperCI", "Beta_Dad_PatGN_P_FDR", "Paternal PGS (Duos in Complete Trios)"),
  pull_series("Paternal_Duo_Full",   "Beta_Dad_PatGN_r_Est", "Beta_Dad_PatGN_r_LowerCI", "Beta_Dad_PatGN_r_UpperCI", "Beta_Dad_PatGN_P_FDR", "Paternal PGS (All Duos)")
) %>%
  mutate(
    Series = factor(Series, levels = p6_levels),
    Stars = fdr_stars(p),
    Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
    Domain       = factor(Domain, levels = domain_order)
  )

dodge_width_p6 <- 0.75

plot6 <- ggplot(p6_data, aes(x = est, y = Unique_Label, fill = Series)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_col(position = position_dodge(width = dodge_width_p6), width = dodge_width_p6 * 0.85) +
  geom_errorbar(aes(xmin = lo, xmax = hi),
                position = position_dodge(width = dodge_width_p6),
                width = dodge_width_p6 * 0.35, linewidth = 0.6, color = "grey25") +
  geom_text(aes(x = hi + max(abs(hi), na.rm = TRUE) * 0.05, label = Stars),
            position = position_dodge(width = dodge_width_p6),
            vjust = 0.75, size = 4.5, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = p6_colors,
                    guide = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE)) +
  scale_y_discrete(
    labels = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  ) +
  scale_x_continuous(breaks = BRK_BETA, expand = expansion(mult = c(0.05, 0.10))) +
  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +
  labs(x = X_BETA, y = NULL, fill = "",
       title = "Sample-composition sensitivity: duo-model genetic nurture estimates in the trio-only versus full duo sample") +
  theme_gn(grid = "x") +
  coord_cartesian(clip = "off")

save_gn(plot6, "RegSup4_duoSubsets", width = 24, height = 24)

p78_est_levels <- c("Child PGS (Population)", "Child PGS (Direct)",
                    "Maternal PGS (Genetic Nurture)", "Paternal PGS (Genetic Nurture)")

build_group_forest <- function(subgroup_top, subgroup_bottom,
                               label_top, label_bottom,
                               colors, plot_title) {

  dat <- read.csv(paste0(outFileStem, "GN_Reg.csv")) %>%
    filter(Subgroup %in% c(subgroup_top, subgroup_bottom))

  grp_levels <- c(label_bottom, label_top)

  dat %>%
    select(Unique_Label, Domain, Subgroup,
           est_pop   = Beta_Child_Pop_Est,
           est_child = Beta_Child_OverallGN_Est,
           est_mat   = Beta_Mom_OverallGN_r_Est,
           est_pat   = Beta_Dad_OverallGN_r_Est,

           lo_pop    = Beta_Child_Pop_LowerCI,
           hi_pop    = Beta_Child_Pop_UpperCI,
           lo_child  = Beta_Child_OverallGN_LowerCI,
           hi_child  = Beta_Child_OverallGN_UpperCI,
           lo_mat    = Beta_Mom_OverallGN_r_LowerCI,
           hi_mat    = Beta_Mom_OverallGN_r_UpperCI,
           lo_pat    = Beta_Dad_OverallGN_r_LowerCI,
           hi_pat    = Beta_Dad_OverallGN_r_UpperCI,

           p_pop     = Beta_Child_Pop_P_FDR,
           p_child   = Beta_Child_OverallGN_P_FDR,
           p_mat     = Beta_Mom_OverallGN_P_FDR,
           p_pat     = Beta_Dad_OverallGN_P_FDR) %>%
    pivot_longer(cols = -c(Unique_Label, Domain, Subgroup),
                 names_to = c(".value", "EstKey"),
                 names_pattern = "(est|lo|hi|p)_(.*)") %>%
    mutate(
      Estimate = dplyr::recode(EstKey,
                               pop   = "Child PGS (Population)",
                               child = "Child PGS (Direct)",
                               mat   = "Maternal PGS (Genetic Nurture)",
                               pat   = "Paternal PGS (Genetic Nurture)"),
      Estimate = factor(Estimate, levels = p78_est_levels),

      Group = ifelse(Subgroup == subgroup_top, label_top, label_bottom),
      Group = factor(Group, levels = grp_levels),

      Stars = fdr_stars(p),

      Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
      Domain       = factor(Domain, levels = domain_order)
    ) -> pdat

  dodge_width_p78 <- 0.62

  ggplot(pdat, aes(x = est, y = Unique_Label, color = Group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_linerange(aes(xmin = lo, xmax = hi),
                   position = position_dodge(width = dodge_width_p78), linewidth = 0.8) +
    geom_point(size = 2.6, position = position_dodge(width = dodge_width_p78)) +
    geom_text(aes(x = hi + max(abs(hi), na.rm = TRUE) * 0.04, label = Stars),
              position = position_dodge(width = dodge_width_p78),
              vjust = 0.75, size = 4.2, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = colors, breaks = grp_levels,
                       guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(
      labels = function(x) stringr::str_wrap(
        paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))), width = 28)
    ) +
    scale_x_continuous(breaks = BRK_BETA, expand = expansion(mult = c(0.05, 0.10))) +
    facet_grid(Domain ~ Estimate, scales = "free_y", space = "free_y") +
    labs(x = X_BETA, y = NULL, color = NULL, title = plot_title) +
    theme_gn(grid = "x") +
    theme(strip.text.y    = element_text(size = GN_BASE_SIZE, face = "bold", angle = -90,
                                         margin = margin(l = 8, r = 8)),
          panel.spacing.y = unit(1.0, "lines"),
          panel.spacing.x = unit(0.8, "lines")) +
    coord_cartesian(clip = "off")
}

calc_height_p78 <- 4.0 + length(trait_order) * 0.55

SEX_COLS <- SEX_COL

plot7 <- build_group_forest(
  subgroup_top    = "Female_Trios", label_top    = "Female Trios",
  subgroup_bottom = "Male_Trios",   label_bottom = "Male Trios",
  colors          = SEX_COLS,
  plot_title      = "Genetic nurture effects by offspring sex"
)

save_gn(plot7, "RegSup5_SexDiff", width = 24, height = calc_height_p78)

ZYG_COLS <- ZYG_COL

plot8 <- build_group_forest(
  subgroup_top    = "MZ_Trios", label_top    = "MZ Trios",
  subgroup_bottom = "DZ_Trios", label_bottom = "DZ Trios",
  colors          = ZYG_COLS,
  plot_title      = "Genetic nurture effects by zygosity"
)

save_gn(plot8, "RegSup6_Zygosity", width = 24, height = calc_height_p78)

forest_colors_sib <- c(
  SERIES4_COL,
  "Sibling PGS (Indirect Genetic Effect)" = SIB_COL
)

GN_Reg_sibIGE <- read.csv(paste0(outFileStem, "GN_Reg.csv")) %>%
  filter(Subgroup == "DZ_SibIGE_Quads")

p9_data <- GN_Reg_sibIGE %>%
  select(Unique_Label, Domain,
         est_pop   = Beta_Child_Pop_Est,
         est_child = Beta_Child_SibGN_Est,
         est_mat   = Beta_Mom_SibGN_r_Est,
         est_pat   = Beta_Dad_SibGN_r_Est,
         est_sib   = Beta_Cotwin_SibGN_Est,

         lo_pop    = Beta_Child_Pop_LowerCI,
         hi_pop    = Beta_Child_Pop_UpperCI,
         lo_child  = Beta_Child_SibGN_LowerCI,
         hi_child  = Beta_Child_SibGN_UpperCI,
         lo_mat    = Beta_Mom_SibGN_r_LowerCI,
         hi_mat    = Beta_Mom_SibGN_r_UpperCI,
         lo_pat    = Beta_Dad_SibGN_r_LowerCI,
         hi_pat    = Beta_Dad_SibGN_r_UpperCI,
         lo_sib    = Beta_Cotwin_SibGN_LowerCI,
         hi_sib    = Beta_Cotwin_SibGN_UpperCI,

         p_pop     = Beta_Child_Pop_P_FDR,
         p_child   = Beta_Child_SibGN_P_FDR,
         p_mat     = Beta_Mom_SibGN_P_FDR,
         p_pat     = Beta_Dad_SibGN_P_FDR,
         p_sib     = Beta_Cotwin_SibGN_P_FDR
  ) %>%
  pivot_longer(
    cols      = -c(Unique_Label, Domain),
    names_to  = c(".value", "Model"),
    names_pattern = "(est|lo|hi|p)_(.*)"
  ) %>%
  mutate(
    Model = dplyr::recode(Model,
                          pop   = "Child PGS (Population)",
                          child = "Child PGS (Direct)",
                          mat   = "Maternal PGS (Genetic Nurture)",
                          pat   = "Paternal PGS (Genetic Nurture)",
                          sib   = "Sibling PGS (Indirect Genetic Effect)"),

    Model = factor(Model, levels = c("Sibling PGS (Indirect Genetic Effect)",
                                     "Paternal PGS (Genetic Nurture)",
                                     "Maternal PGS (Genetic Nurture)",
                                     "Child PGS (Direct)",
                                     "Child PGS (Population)")),

    trait_clean = paste0(toupper(substr(Unique_Label, 1, 1)), substr(Unique_Label, 2, nchar(Unique_Label))),

    Unique_Label = factor(Unique_Label, levels = rev(trait_order)),
    trait_clean  = factor(trait_clean,
                          levels = unique(trait_clean[order(match(Unique_Label, levels(Unique_Label)))])),
    Domain       = factor(Domain, levels = domain_order),

    Stars = fdr_stars(p)
  )

dodge_width_p9 <- 0.7

plot9 <- ggplot(p9_data, aes(x = est, y = trait_clean, color = Model)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey50") +

  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = dodge_width_p9), linewidth = 0.8) +

  geom_point(size = 3, position = position_dodge(width = dodge_width_p9)) +

  geom_text(aes(x = hi + max(abs(hi), na.rm = TRUE) * 0.05, label = Stars),
            position = position_dodge(width = dodge_width_p9),
            vjust = 0.75, size = 4.5, fontface = "bold", show.legend = FALSE) +

  scale_color_manual(values = forest_colors_sib,
                     breaks = c("Child PGS (Population)", "Child PGS (Direct)",
                                "Maternal PGS (Genetic Nurture)", "Paternal PGS (Genetic Nurture)",
                                "Sibling PGS (Indirect Genetic Effect)")) +
  scale_x_continuous(breaks = BRK_BETA) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +

  facet_wrap(~ Domain, scales = "free_y", ncol = 2) +

  labs(
    title    = "Child, parental, and sibling PGS effects (DZ co-twin sensitivity)",
    x        = X_BETA,
    y        = NULL,
    color    = NULL
  ) +
  theme_gn(grid = "x") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

calc_height_p9 <- 4.0 + ((length(unique(p9_data$Unique_Label)) / 2) * 0.8)

save_gn(plot9, "RegSup7_SiblingIGE", width = 24, height = calc_height_p9)

