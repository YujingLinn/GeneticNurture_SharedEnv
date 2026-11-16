library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)

regFileStem <- "GN_Results/Reg/"
aceFileStem <- "GN_Results/TwinResults/"
outFileStem <- "GN_Results/ACE_vs_GN/"
if (!dir.exists(outFileStem)) dir.create(outFileStem, recursive = TRUE)

source("0GN_Theme.R")
MOM_LAB <- "Maternal Genetic Nurture"
DAD_LAB <- "Paternal Genetic Nurture"
PAR_COL <- setNames(c(MAT_COL, PAT_COL), c(MOM_LAB, DAD_LAB))

df_reg <- read.csv(paste0(regFileStem, "GN_Reg.csv"))
ace    <- read_excel(paste0(aceFileStem, "uniGN_121026.xlsx"), sheet = "ParentDNA")

ace_t <- data.frame(trait = ace$Specific_Label, Domain_ace = ace$Category,
                    c2 = ace$c2,
                    stringsAsFactors = FALSE) %>%
  mutate(c = sqrt(pmax(c2, 0)))

reg <- df_reg %>% filter(Subgroup == "All_Trios")

reg_t <- reg %>%
  transmute(trait  = Unique_Label, Domain,
            pop    = Beta_Child_Pop_Est,
            mom    = Beta_Mom_OverallGN_r_Est,
            dad    = Beta_Dad_OverallGN_r_Est)

tr <- inner_join(reg_t, ace_t, by = "trait")

dropped <- tr$trait[is.na(tr$c2) | is.na(tr$mom) | is.na(tr$dad) | is.na(tr$pop)]
if (length(dropped))
  message(sprintf("Dropped %d traits (ACE not converged or GN/Pop missing): %s",
                  length(dropped), paste(dropped, collapse = ", ")))
tr <- tr %>% filter(!trait %in% dropped) %>% select(-Domain_ace)

GROUP_LEVELS <- c("g", "Verbal Ability", "Nonverbal Ability",
                  "Educational Achievement", "Educational Attainment",
                  "Height", "BMI", "Psychopathology")

tr <- tr %>%
  mutate(construct = trimws(sub("\\s+age\\s+[0-9]+$", "", trait)),
         age       = suppressWarnings(as.numeric(sub(".*age\\s+([0-9]+)$", "\\1", trait))),
         group = case_when(
           Domain == "Psychopathology"                                ~ "Psychopathology",
           construct == "g"                                           ~ "g",
           construct == "verbal ability"                              ~ "Verbal Ability",
           construct == "nonverbal ability"                           ~ "Nonverbal Ability",
           construct %in% c("primary school grades", "GCSE grades",
                            "A-level grades", "university grades")    ~ "Educational Achievement",
           construct %in% c("A-level enrolment", "university enrolment",
                            "years of schooling")                     ~ "Educational Attainment",
           construct == "height"                                      ~ "Height",
           construct == "BMI"                                         ~ "BMI",
           TRUE                                                       ~ NA_character_))
tr$group <- factor(tr$group, levels = GROUP_LEVELS)

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

pd_grp <- tr %>%
  left_join(trait_pgs, by = "trait") %>%
  left_join(fam_cor, by = "pgs_label") %>%
  mutate(gn_pop = r_mc * mom + r_fc * dad) %>%
  group_by(group) %>%
  summarise(n_traits    = n(),
            c_mean      = mean(c),
            gn_pop_mean = mean(gn_pop),
            .groups     = "drop")

st6 <- pd_grp %>%
  summarise(r = cor(c_mean, gn_pop_mean), p = cor.test(c_mean, gn_pop_mean)$p.value, n = n()) %>%
  mutate(lab = sprintf("  r = %.2f (p = %.3f, n = %d)", r, p, n))

p6 <- ggplot(pd_grp, aes(x = c_mean, y = gn_pop_mean)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9,
              colour = OVR_COL) +
  geom_point(size = 3.8, colour = OVR_COL) +
  ggrepel::geom_text_repel(aes(label = group), size = 4.2, colour = "grey30",
                           box.padding = 0.4, point.padding = 0.3,
                           min.segment.length = 0.3, segment.colour = "grey70", seed = 1) +
  geom_text(data = st6, aes(x = -Inf, y = Inf, label = lab), vjust = 1.5, inherit.aes = FALSE,
            hjust = 0, size = 5.5, colour = "grey25") +
  scale_x_continuous(breaks = BRK_BETA, expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.02),
                     expand = expansion(mult = c(0.08, 0.22))) +
  labs(x = expression("Shared environment: " * italic(c) ~ "path coefficient"),
       y = "Genetic nurture",
       title = "Genetic nurture and shared environment across trait groups") +
  theme_gn(grid = "both")
save_gn(p6, "ACE_vs_GN_Domain_PopMinusDirect", width = 9.5, height = 7.5)

RATIO_COL <- "#B5446E"

rt_grp <- tr %>%
  left_join(trait_pgs, by = "trait") %>%
  left_join(fam_cor, by = "pgs_label") %>%
  mutate(gn_pop = r_mc * mom + r_fc * dad) %>%
  group_by(group) %>%
  summarise(n_traits   = n(),
            c_mean     = mean(c),
            ratio_mean = mean(gn_pop) / mean(pop),
            .groups    = "drop")

st7 <- rt_grp %>%
  summarise(r = cor(c_mean, ratio_mean), p = cor.test(c_mean, ratio_mean)$p.value, n = n()) %>%
  mutate(lab = sprintf("  r = %.2f (p = %.3f, n = %d)", r, p, n))

p7 <- ggplot(rt_grp, aes(x = c_mean, y = ratio_mean)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9,
              colour = RATIO_COL) +
  geom_point(size = 3.8, colour = RATIO_COL) +
  ggrepel::geom_text_repel(aes(label = group), size = 4.2, colour = "grey30",
                           box.padding = 0.4, point.padding = 0.3,
                           min.segment.length = 0.3, segment.colour = "grey70", seed = 1) +
  geom_text(data = st7, aes(x = -Inf, y = Inf, label = lab), vjust = 1.5, inherit.aes = FALSE,
            hjust = 0, size = 5.5, colour = "grey25") +
  scale_x_continuous(breaks = BRK_BETA, expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(breaks = seq(-1, 2, by = 0.05),
                     expand = expansion(mult = c(0.08, 0.22))) +
  labs(x = expression("Shared environment: " * italic(c) ~ "path coefficient"),
       y = "Genetic nurture / population prediction",
       title = "Genetic nurture and shared environment across trait groups") +
  theme_gn(grid = "both")
save_gn(p7, "ACE_vs_GN_Domain_GNPopRatio", width = 9.5, height = 7.5)

par_grp <- tr %>%
  left_join(trait_pgs, by = "trait") %>%
  left_join(fam_cor, by = "pgs_label") %>%
  mutate(gn_m = r_mc * mom,
         gn_f = r_fc * dad) %>%
  group_by(group) %>%
  summarise(n_traits = n(),
            c_mean   = mean(c),
            pop_mean = mean(pop),
            gn_m     = mean(gn_m),
            gn_f     = mean(gn_f),
            .groups  = "drop")

par_long <- bind_rows(
  par_grp %>% transmute(group, c_mean, parent = MOM_LAB, gn_pop = gn_m, ratio = gn_m / pop_mean),
  par_grp %>% transmute(group, c_mean, parent = DAD_LAB, gn_pop = gn_f, ratio = gn_f / pop_mean)) %>%
  mutate(parent = factor(parent, levels = c(MOM_LAB, DAD_LAB)))

plot_par_xy <- function(d, yvar, ylab, ybreaks, fname) {
  d$y <- d[[yvar]]
  st  <- d %>% group_by(parent) %>%
    summarise(r = cor(c_mean, y), p = cor.test(c_mean, y)$p.value, n = n(), .groups = "drop") %>%
    mutate(lab = sprintf("  r = %.2f (p = %.3f, n = %d)", r, p, n))
  p <- ggplot(d, aes(x = c_mean, y = y, colour = parent)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9,
                show.legend = FALSE) +
    geom_point(size = 3.8, show.legend = FALSE) +
    ggrepel::geom_text_repel(aes(label = group), size = 4.2, colour = "grey30",
                             box.padding = 0.4, point.padding = 0.3,
                             min.segment.length = 0.3, segment.colour = "grey70", seed = 1) +
    geom_text(data = st, aes(x = -Inf, y = Inf, label = lab), vjust = 1.5, inherit.aes = FALSE,
              hjust = 0, size = 5.5, colour = "grey25") +
    facet_wrap(~ parent, nrow = 1) +
    scale_colour_manual(values = PAR_COL, guide = "none") +
    scale_x_continuous(breaks = BRK_BETA, expand = expansion(mult = c(0.08, 0.08))) +
    scale_y_continuous(breaks = ybreaks, expand = expansion(mult = c(0.08, 0.22))) +
    labs(x = expression("Shared environment: " * italic(c) ~ "path coefficient"), y = ylab,
         title = "Genetic nurture and shared environment across trait groups") +
    theme_gn(grid = "both")
  save_gn(p, fname, width = 15, height = 7.5)
  invisible(p)
}

p9  <- plot_par_xy(par_long, "gn_pop", "Genetic nurture",
                   seq(-1, 1, by = 0.01), "ACE_vs_GN_Domain_PopMinusDirect_ByParent")
p10 <- plot_par_xy(par_long, "ratio", "Genetic nurture / population prediction",
                   seq(-1, 2, by = 0.05), "ACE_vs_GN_Domain_GNPopRatio_ByParent")
