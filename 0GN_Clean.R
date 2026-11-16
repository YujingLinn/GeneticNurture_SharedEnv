library(plyr)
library(tidyverse)
library(purrr)
library(Hmisc)
library(foreign)
library(stringr)
library(readr)
library(dplyr)

library(conflicted)
conflicts_prefer(
  dplyr::select, dplyr::mutate, dplyr::filter, dplyr::group_by, dplyr::ungroup,
  dplyr::summarise, dplyr::summarize, dplyr::rename, dplyr::arrange, dplyr::distinct,
  dplyr::count, dplyr::first, dplyr::last, dplyr::desc,
  dplyr::left_join, dplyr::right_join, dplyr::full_join, dplyr::inner_join
)

GN_base     <- "."
spss_file   <- file.path(GN_base, "GN_Data", "file_name.sav")

PGS_dir     <- file.path(GN_base, "PGS")
pca_file    <- file.path(GN_base, "teds_trio_projected_pca.sscore.txt")
outFileStem <- file.path(GN_base, "GN_Data")

all_files <- list.files(
  path = PGS_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)
process_pgs_file <- function(file_path) {

  pgs_name <- tools::file_path_sans_ext(basename(file_path))

  score_column <- if (str_detect(pgs_name, "pred_auto")) {
    "final_pred_auto"
  } else if (str_detect(pgs_name, "pred_inf")) {
    "pred_inf"
  } else {
    stop(paste("Unknown score type in filename:", pgs_name))
  }

  child_score_col  <- pgs_name
  mother_score_col <- paste0(pgs_name, "_fp")
  father_score_col <- paste0(pgs_name, "_mp")

  df <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(
      family.ID   = as.character(family.ID),
      sample.ID   = as.character(sample.ID),
      member_code = str_remove(sample.ID, paste0("^", family.ID))
    )

  children <- df %>%
    filter(member_code %in% c("1", "2")) %>%
    mutate(id_fam = as.numeric(family.ID),
           id_twin = as.numeric(sample.ID)) %>%
    select(id_fam, id_twin, !!child_score_col := all_of(score_column))

  mothers_scores <- df %>%
    filter(member_code == "3") %>%
    distinct(family.ID, .keep_all = TRUE) %>%
    mutate(id_fam = as.numeric(family.ID)) %>%
    select(id_fam, !!mother_score_col := all_of(score_column))

  fathers_scores <- df %>%
    filter(member_code == "4") %>%
    distinct(family.ID, .keep_all = TRUE) %>%
    mutate(id_fam = as.numeric(family.ID)) %>%
    select(id_fam, !!father_score_col := all_of(score_column))

  children %>%
    left_join(mothers_scores, by = "id_fam") %>%
    left_join(fathers_scores, by = "id_fam")
}

list_of_processed_dfs <- map(all_files, process_pgs_file)

pgs_data <- list_of_processed_dfs %>%
  reduce(full_join, by = c("id_fam", "id_twin"))

pca_projected <- read.table(pca_file, header = TRUE, comment.char = "")
colnames(pca_projected)[1] <- "id_fam"
colnames(pca_projected)[2] <- "id_twin"

projected_scores <- pca_projected[, c("id_fam", "id_twin", paste0("SCORE", 1:20, "_AVG"))] %>%
  rename_with(~ paste0("PC", seq_along(.x)), .cols = starts_with("SCORE")) %>%
  mutate(
    id_fam      = as.character(id_fam),
    id_twin     = as.character(id_twin),
    member_code = str_remove(id_twin, paste0("^", id_fam))
  )

pcs_twins <- projected_scores %>%
  filter(member_code %in% c("1", "2")) %>%
  mutate(id_fam = as.numeric(id_fam), id_twin = as.numeric(id_twin)) %>%
  select(-member_code)

pcs_mothers <- projected_scores %>%
  filter(member_code == "3") %>%
  distinct(id_fam, .keep_all = TRUE) %>%
  mutate(id_fam = as.numeric(id_fam)) %>%
  rename_with(~ paste0(.x, "_fp"), .cols = starts_with("PC")) %>%
  select(-id_twin, -member_code)

pcs_fathers <- projected_scores %>%
  filter(member_code == "4") %>%
  distinct(id_fam, .keep_all = TRUE) %>%
  mutate(id_fam = as.numeric(id_fam)) %>%
  rename_with(~ paste0(.x, "_mp"), .cols = starts_with("PC")) %>%
  select(-id_twin, -member_code)

updated_PCA <- pcs_twins %>%
  left_join(pcs_mothers, by = "id_fam") %>%
  left_join(pcs_fathers, by = "id_fam")

updated_data <- pgs_data %>%
  left_join(updated_PCA, by = c("id_fam", "id_twin"))

pc_child  <- paste0("PC", 1:10)
pc_mother <- paste0("PC", 1:10, "_fp")
pc_father <- paste0("PC", 1:10, "_mp")

resid_child_cols <- grep("_pred_(auto|inf)$", names(updated_data), value = TRUE)

for (cc in resid_child_cols) {
  for (col in c(cc, paste0(cc, "_fp"), paste0(cc, "_mp"))) {
    if (col %in% names(updated_data))
      updated_data[[paste0(col, "_noPC")]] <- updated_data[[col]]
  }
}

.resid_pc_self <- function(df, y_col, pc_cols) {
  tmp <- data.frame(y = df[[y_col]], df[, pc_cols, drop = FALSE])
  fit <- lm(y ~ ., data = tmp, na.action = na.exclude)
  as.numeric(residuals(fit))
}
.resid_pc_family <- function(df, y_col, pc_cols, fam_col = "id_fam") {
  d1  <- df[!duplicated(df[[fam_col]]), c(y_col, pc_cols)]
  tmp <- data.frame(y = d1[[y_col]], d1[, pc_cols, drop = FALSE])
  fit <- lm(y ~ ., data = tmp)
  nd  <- df[, pc_cols, drop = FALSE]
  colnames(nd) <- colnames(d1)[-1]
  pred <- as.numeric(predict(fit, newdata = nd))
  df[[y_col]] - pred
}

for (cc in resid_child_cols) {
  fp <- paste0(cc, "_fp")
  mp <- paste0(cc, "_mp")

  updated_data[[cc]] <- .resid_pc_self(updated_data, cc, pc_child)
  if (fp %in% names(updated_data))
    updated_data[[fp]] <- .resid_pc_family(updated_data, fp, pc_mother)
  if (mp %in% names(updated_data))
    updated_data[[mp]] <- .resid_pc_family(updated_data, mp, pc_father)
}

origin_LDPred1_IQ3 <- read.spss(spss_file, to.data.frame = TRUE, use.value.labels = FALSE)
matching_cols <- origin_LDPred1_IQ3 %>%
  select(starts_with(c(paste0("PC", 1:10),
                       "EA4_no23andme_Okbay2022", "IQ_Savage2018_fp", "IQ_Savage2018_mp",
                       "ChildhoodBMI", "BMI_Giant", "Height_Yengo",
                       "Depression", "Anxiety", "ADHD", "ASD_Grove", "AnorexiaNervosa",
                       "BipolarDisorder", "Stress", "Alcohol", "ExecutiveFunction",
                       "Cognitive_performance_Lee", "Cognitive_Demange"))) %>%
  names()

origin_LDPred1_IQ3_reduced <- origin_LDPred1_IQ3 %>%
  select(-all_of(matching_cols))

origin_LDPred1_IQ3_reduced_updated <- full_join(
  origin_LDPred1_IQ3_reduced,
  updated_data,
  by = c("id_fam", "id_twin")
)

GN_raw <- origin_LDPred1_IQ3_reduced_updated[origin_LDPred1_IQ3_reduced_updated$exclude1 == 0, ]

df_verbal12T1    <- select(GN_raw, id_twin, ltotota1, lfltot1, lintot1, lpitota1, lgotota1, lyntot1, ltowt1)
df_verbal12T2    <- select(GN_raw, id_twin, ltotota2, lfltot2, lintot2, lpitota2, lgotota2, lyntot2, ltowt2)
df_nonverbal12T1 <- select(GN_raw, id_twin, lma1tot1, lma2tot1, lma3tot1)
df_nonverbal12T2 <- select(GN_raw, id_twin, lma1tot2, lma2tot2, lma3tot2)

df_list_cognitive12 <- list(df_verbal12T1, df_verbal12T2, df_nonverbal12T1, df_nonverbal12T2)

score_list_cognitive12 <- map(df_list_cognitive12, ~ {
  n_missing <- apply(.[-1], 1, function(x) sum(is.na(x))) > ((ncol(.) - 1) / 2)
  row_sums  <- apply(.[!n_missing, -1], 1, mean, na.rm = TRUE) * (ncol(.) - 1)
  data.frame(
    id_twin  = .[!n_missing, 1],
    row_sums = row_sums
  )
})

df_cognitive12 <- score_list_cognitive12 %>% reduce(full_join, by = "id_twin")
colnames(df_cognitive12) <- c("id_twin", "lverbal12T1", "lverbal12T2", "lnonverbal12T1", "lnonverbal12T2")

GN_raw <- GN_raw %>% full_join(df_cognitive12)
NumberGradesReport_a_level_enrolment_conversion <- function(data, a_level_reported, as_level_reported, new_item) {
  data[[new_item]] <- ifelse(
    !is.na(data[[a_level_reported]]),
    ifelse(data[[a_level_reported]] %in% 1:6, 1,
           ifelse(data[[a_level_reported]] == 0,
                  ifelse(data[[as_level_reported]] %in% 1:5, 1,
                         ifelse(data[[as_level_reported]] == 0, 0, NA)),
                  NA)),
    ifelse(!is.na(data[[as_level_reported]]),
           ifelse(data[[as_level_reported]] %in% 1:5, 1,
                  ifelse(data[[as_level_reported]] == 0, 0, NA)),
           NA)
  )
  return(data)
}

GN_raw <- NumberGradesReport_a_level_enrolment_conversion(GN_raw, "rcqaln1", "rcqasn1", "ra_level_enrol1")
GN_raw <- NumberGradesReport_a_level_enrolment_conversion(GN_raw, "rcqaln2", "rcqasn2", "ra_level_enrol2")

edu_attain_collapse <- function(data, EA26, EA21, new_item) {
  data[[new_item]] <- ifelse(
    !is.na(data[[EA26]]),
    case_when(
      data[[EA26]] %in% 1      ~ NA,
      data[[EA26]] %in% c(2,3,4) ~ 10,
      data[[EA26]] %in% c(5,6) ~ 12,
      data[[EA26]] %in% c(7,8) ~ 14,
      data[[EA26]] == 9        ~ 15,
      data[[EA26]] == 10       ~ 16,
      data[[EA26]] == 11       ~ 20
    ),
    ifelse(
      !is.na(data[[EA21]]),
      case_when(
        data[[EA21]] %in% 1      ~ NA,
        data[[EA21]] %in% c(2,3,4) ~ 10,
        data[[EA21]] %in% c(5,6) ~ 12,
        data[[EA21]] %in% c(7,8) ~ 14,
        data[[EA21]] == 9        ~ 15,
        data[[EA21]] == 10       ~ 16,
        data[[EA21]] == 11       ~ 20
      ),
      NA
    )
  )
  return(data)
}

GN_raw <- edu_attain_collapse(GN_raw, "zmhhqual1", "u1chqualp1", "zEA1")
GN_raw <- edu_attain_collapse(GN_raw, "zmhhqual2", "u1chqualp2", "zEA2")

GN_raw$zEA1 <- as.numeric(GN_raw$zEA1)
GN_raw$zEA2 <- as.numeric(GN_raw$zEA2)

df <- GN_raw

calculate_prorated_average <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  total_items     <- ncol(x)
  available_items <- rowSums(!is.na(x))
  average_score   <- rowMeans(x, na.rm = TRUE)
  ifelse(available_items >= (total_items / 2), average_score, NA)
}

cols_eat1 <- c("zmheatd011",
               "zmheatd03a1", "zmheatd03b1", "zmheatd03c1", "zmheatd03d1",
               "zmheatd041",
               "zmheatd05a1", "zmheatd05b1", "zmheatd05c1",
               "zmheatd06a1", "zmheatd06b1",
               "zmheatd06c1", "zmheatd06d1", "zmheatd06e1",
               "zmheatd071",
               "zmheatd09a1", "zmheatd09b1", "zmheatd09c1", "zmheatd09d1",
               "zmheatd101",
               "zmheatd11a1", "zmheatd11b1", "zmheatd11c1", "zmheatd11d1")

cols_eat2 <- c("zmheatd012",
               "zmheatd03a2", "zmheatd03b2", "zmheatd03c2", "zmheatd03d2",
               "zmheatd042",
               "zmheatd05a2", "zmheatd05b2", "zmheatd05c2",
               "zmheatd06a2", "zmheatd06b2",
               "zmheatd06c2", "zmheatd06d2", "zmheatd06e2",
               "zmheatd072",
               "zmheatd09a2", "zmheatd09b2", "zmheatd09c2", "zmheatd09d2",
               "zmheatd102",
               "zmheatd11a2", "zmheatd11b2", "zmheatd11c2", "zmheatd11d2")

after_01 <- c("03a", "03b", "03c", "03d")
after_04 <- c("05a", "05b", "05c", "06a", "06b", "06c", "06d", "06e", "07",
              "09a", "09b", "09c", "09d")

for (s in 1:2) {
  scr01 <- !is.na(df[[paste0("zmheatd01", s)]]) & df[[paste0("zmheatd01", s)]] == 0
  scr04 <- !is.na(df[[paste0("zmheatd04", s)]]) & df[[paste0("zmheatd04", s)]] == 0
  for (v in paste0("zmheatd", after_01, s)) df[[v]][scr01 & is.na(df[[v]])] <- 0
  for (v in paste0("zmheatd", after_04, s)) df[[v]][scr04 & is.na(df[[v]])] <- 0
}

df$zmheatdt1 <- calculate_prorated_average(df[, cols_eat1])
df$zmheatdt2 <- calculate_prorated_average(df[, cols_eat2])

GN_raw <- GN_raw %>%
  left_join(df %>% select(id_twin, id_fam, zmheatdt1, zmheatdt2), by = c("id_twin", "id_fam"))

GN_raw$zmhmdq_bipolar1 <- GN_raw$zmhmdq1count1
GN_raw$zmhmdq_bipolar2 <- GN_raw$zmhmdq1count2

GN_raw_full   <- GN_raw
ref_unrelated <- GN_raw_full %>% filter(selectunpaired == 1)

GN_scaled_full <- GN_raw_full

child_pgs_cols <- grep("_pred_(auto|inf)$", names(GN_raw_full), value = TRUE)

for (child_col in child_pgs_cols) {
  fp_col <- paste0(child_col, "_fp")
  mp_col <- paste0(child_col, "_mp")

  mu      <- mean(ref_unrelated[[child_col]], na.rm = TRUE)
  sd_c    <- sd(ref_unrelated[[child_col]],   na.rm = TRUE)

  GN_scaled_full[[child_col]] <- (GN_scaled_full[[child_col]] - mu) / sd_c
  if (fp_col %in% names(GN_scaled_full))
    GN_scaled_full[[fp_col]] <- (GN_scaled_full[[fp_col]] - mu) / sd_c
  if (mp_col %in% names(GN_scaled_full))
    GN_scaled_full[[mp_col]] <- (GN_scaled_full[[mp_col]] - mu) / sd_c
}

for (child_col in child_pgs_cols) {
  npc_c  <- paste0(child_col, "_noPC")
  npc_fp <- paste0(child_col, "_fp_noPC")
  npc_mp <- paste0(child_col, "_mp_noPC")

  mu   <- mean(ref_unrelated[[npc_c]], na.rm = TRUE)
  sd_c <- sd(ref_unrelated[[npc_c]],   na.rm = TRUE)

  GN_scaled_full[[npc_c]] <- (GN_scaled_full[[npc_c]] - mu) / sd_c
  if (npc_fp %in% names(GN_scaled_full))
    GN_scaled_full[[npc_fp]] <- (GN_scaled_full[[npc_fp]] - mu) / sd_c
  if (npc_mp %in% names(GN_scaled_full))
    GN_scaled_full[[npc_mp]] <- (GN_scaled_full[[npc_mp]] - mu) / sd_c
}

outcome_vars <- c(
  "gcg","icg","jcg","lcg","pcg","ucgt",
  "gcl","icvb","jcvb","lverbal12T","pcvctota","ucgvbt",
  "gcn","icnv","jcnv","lnonverbal12T","pcrvtota","ucgnvt",
  "gt2ac","it3ac","jt3ac","lt3ac","pcexgcsecoregrdm","rcqalgrdm","u1cdegr1",
  "ra_level_enrol","rcqhe","zEA",
  "ghtcm","lchtcm","nchtcm","pcqdhtcm","u1chtcm","zmhheight",
  "gbmi","lbmi","ncbmi","pcbmi","u1cbmi","zmhbmi",
  "zmhpclt","zmhmfqt","zmhalcoaudit","zmhganxt","zmheatdt","zmhconnt","zmhmdq_bipolar","zmhraadst"
)

pheno_vars <- c(paste0(outcome_vars, "1"), paste0(outcome_vars, "2"))

vars_to_scale <- pheno_vars

vars_to_scale <- intersect(vars_to_scale, names(GN_raw_full))

unrel_stats <- ref_unrelated %>%
  summarise(across(all_of(vars_to_scale),
                   list(mean = ~mean(., na.rm = TRUE),
                        sd   = ~sd(.,   na.rm = TRUE))))

GN_scaled_full <- GN_scaled_full %>%
  mutate(across(all_of(vars_to_scale),
                ~ (. - unrel_stats[[paste0(cur_column(), "_mean")]]) /
                  unrel_stats[[paste0(cur_column(), "_sd")]]))

GN_raw_selectunpaired <- subset(GN_raw_full,    selectunpaired == 1)

GN_scaled             <- GN_scaled_full
GN_scaled_selectunpaired <- subset(GN_scaled_full, selectunpaired == 1)

write.csv(GN_raw_full,           file.path(outFileStem, "GN_raw.csv"),                row.names = FALSE)
write.csv(GN_raw_selectunpaired, file.path(outFileStem, "GN_raw_selectunpaired.csv"), row.names = FALSE)

write.csv(GN_scaled,                file.path(outFileStem, "GN_scaled.csv"),                row.names = FALSE)
write.csv(GN_scaled_selectunpaired, file.path(outFileStem, "GN_scaled_selectunpaired.csv"), row.names = FALSE)
