library(dplyr)
library(openxlsx)

outFileStem <- "GN_Results/Reg/"

GN_Reg <- read.csv(paste0(outFileStem, "GN_Reg.csv")) %>%
  filter(Subgroup == "All_Trios")

fam_cor <- read.csv("GN_Results/Data_Checking/GN_PGS_family_correlations_trios.csv") %>%
  filter(Trait != "N (pairs)") %>%
  transmute(pgs_label = Trait,
            r_mc = r_MotherChild,
            r_fc = r_FatherChild,
            r_mf = r_MotherFather_AM)

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
construct_order <- c("Educational attainment", "Educational achievement",
                     "General cognitive ability (g)", "Verbal ability", "Nonverbal ability",
                     "Psychopathology", "Height", "BMI")

tab_trait <- GN_Reg %>%
  transmute(trait  = Unique_Label,
            Domain,
            N_used, N_fam,
            pop    = Beta_Child_Pop_Est,
            direct = Beta_Child_OverallGN_Est,
            mom    = Beta_Mom_OverallGN_r_Est,
            dad    = Beta_Dad_OverallGN_r_Est,
            R2_Child_Pop_Est, R2_Covar_Only_Est) %>%
  left_join(trait_pgs, by = "trait") %>%
  left_join(fam_cor, by = "pgs_label") %>%
  mutate(construct = ifelse(Domain == "Psychopathology", "Psychopathology",
                            unname(construct_map[trait])),
         mom_w         = r_mc * mom,
         dad_w         = r_fc * dad,
         direct_share  = direct / pop,
         gn_share      = (mom_w + dad_w) / pop,
         R2_pop        = R2_Child_Pop_Est - R2_Covar_Only_Est,
         V_direct      = direct^2,
         V_GN          = mom^2 + dad^2,
         V_GN_pop      = mom_w^2 + dad_w^2) %>%
  select(construct, trait, Domain, pgs_label, N_used, N_fam,
         r_mc, r_fc, pop, direct, mom, dad,
         mom_w, dad_w, direct_share, gn_share,
         R2_Child_Pop_Est, R2_Covar_Only_Est, R2_pop,
         V_direct, V_GN, V_GN_pop)

tab_group <- tab_trait %>%
  group_by(construct) %>%
  summarise(n_traits = n(),
            across(c(pop, direct, mom, dad, mom_w, dad_w, R2_pop), mean, .names = "mean_{.col}"),
            .groups = "drop") %>%
  mutate(direct_share  = mean_direct / mean_pop,
         gn_share      = (mean_mom_w + mean_dad_w) / mean_pop,
         V_direct      = mean_direct^2,
         V_GN          = mean_mom^2 + mean_dad^2,
         V_GN_pop      = mean_mom_w^2 + mean_dad_w^2) %>%
  arrange(match(construct, construct_order))

s3_group <- tab_group %>%
  transmute(`Trait group`               = construct,
            Outcome                     = paste0("Domain mean (", n_traits, " outcomes)"),
            `Population β`              = mean_pop,
            `Direct β`                  = mean_direct,
            `Maternal GN β`             = mean_mom,
            `Paternal GN β`             = mean_dad,
            `Maternal GN β (weighted)`  = mean_mom_w,
            `Paternal GN β (weighted)`  = mean_dad_w,
            `Direct share`              = direct_share,
            `GN share`                  = gn_share,
            `Population ΔR²`            = mean_R2_pop,
            `Variance: direct`          = V_direct,
            `Variance: GN (total)`      = V_GN,
            `Variance: GN (population)` = V_GN_pop)
s3_trait <- tab_trait %>%
  transmute(`Trait group`               = construct,
            Outcome                     = trait,
            `Population β`              = pop,
            `Direct β`                  = direct,
            `Maternal GN β`             = mom,
            `Paternal GN β`             = dad,
            `Maternal GN β (weighted)`  = mom_w,
            `Paternal GN β (weighted)`  = dad_w,
            `Direct share`              = direct_share,
            `GN share`                  = gn_share,
            `Population ΔR²`            = R2_pop,
            `Variance: direct`          = V_direct,
            `Variance: GN (total)`      = V_GN,
            `Variance: GN (population)` = V_GN_pop)
tab_s3 <- bind_rows(s3_group, s3_trait)

s3_title <- "Supplementary Table 3. Decomposition of population polygenic score prediction into direct and genetic nurture effects"

wb <- createWorkbook()
addWorksheet(wb, "S3")
writeData(wb, "S3", s3_title, startRow = 1, startCol = 1)
addStyle(wb, "S3", createStyle(fontSize = 12, textDecoration = "bold"), rows = 1, cols = 1)
writeData(wb, "S3", tab_s3, startRow = 3, startCol = 1,
          headerStyle = createStyle(textDecoration = "bold", wrapText = TRUE,
                                    border = "bottom", borderStyle = "thin"))
addStyle(wb, "S3", createStyle(numFmt = "0.000"),
         rows = 4:(3 + nrow(tab_s3)), cols = 3:ncol(tab_s3), gridExpand = TRUE)
writeData(wb, "S3", data.frame(note = s3_notes), startRow = 3 + nrow(tab_s3) + 3,
          startCol = 1, colNames = FALSE)
setRowHeights(wb, "S3", rows = 3, heights = 45)
setColWidths(wb, "S3", cols = 1, widths = 40)
setColWidths(wb, "S3", cols = 2:ncol(tab_s3), widths = 14)
freezePane(wb, "S3", firstActiveRow = 4, firstActiveCol = 2)
saveWorkbook(wb, paste0(outFileStem, "GN_Reg_VarianceTables.xlsx"), overwrite = TRUE)
