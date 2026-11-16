library(OpenMx)
library(Hmisc)
library(dplyr)

sourceFileStem <- "GN_Data/"
outFileStem    <- "GN_Results/TwinResults/"

source("./ISG Workshop/ISGW2022/day2/powerFun.R") # post-hoc power analyses

dat_full <- read.csv(paste0(sourceFileStem, "GN_raw.csv"))
dat_whole <- subset(dat_full, random == 0)
pc_child  <- "PC1"
pc_mother <- "PC1_fp"
pc_father <- "PC1_mp"

.needed_pc <- c(pc_child, pc_mother, pc_father)
if (!all(.needed_pc %in% names(dat_whole))) {
  stop("Parent-DNA filter columns not found: ",
       paste(setdiff(.needed_pc, names(dat_whole)), collapse = ", "))
}

keep_parentDNA <- !is.na(dat_whole[[pc_child]]) &
  (!is.na(dat_whole[[pc_mother]]) | !is.na(dat_whole[[pc_father]]))
dat_parentDNA <- dat_whole[keep_parentDNA, ]
cat("Parent-DNA subsample: kept", nrow(dat_parentDNA), "of", nrow(dat_whole), "rows\n")

sample_list <- list(whole = dat_whole, parentDNA = dat_parentDNA)
sample_dirs <- c(whole = "WholeSample/", parentDNA = "ParentDNA/")
for (.d in sample_dirs) dir.create(paste0(outFileStem, .d), recursive = TRUE, showWarnings = FALSE)

min_pairs <- 150

dfBin <- rbind(
  c("ra_level_enrol1", "ra_level_enrol2", "zygos", "rcqalage1", "rcqalage2", "sex1", "sex2"),
  c("rcqhe1", "rcqhe2", "zygos", "rcqalage1", "rcqalage2", "sex1", "sex2")
)

uniGN_all_list <- list()

# ----------------------------------------------------------------------------------------------------------------------
#  Script: 00_ACEvc.R  
#  Author: Sarah Medland
#    Date: 29 02 2020 
#  
# Twin Univariate ACE model to estimate causes of variation across multiple groups
# Matrix style model 
# https://www.colorado.edu/ibg/isg-workshop/online-lectures-and-practicals/e2-practical-ace-models
# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|

for (sample_tag in names(sample_list)) {

  dat    <- sample_list[[sample_tag]]
  outDir <- paste0(outFileStem, sample_dirs[[sample_tag]])
  cat("\n\n########## SAMPLE:", sample_tag, "(N rows =", nrow(dat), ") ##########\n\n")

  uniGN_binary <- data.frame(Variable = character(), Label = character(),
                             nMZ_pairs = numeric(), nDZ_pairs = numeric(),
                             rMZ = numeric(), rDZ = numeric(),
                             h2=numeric(), h2.CI.lower=numeric(), h2.CI.upper=numeric(),
                             c2=numeric(), c2.CI.lower=numeric(), c2.CI.upper=numeric(),
                             e2=numeric(), e2.CI.lower=numeric(), e2.CI.upper=numeric(),
                             stringsAsFactors = FALSE)

  for (i in 1:nrow(dfBin)) {
    selVars <- c(dfBin[i,1], dfBin[i,2])

    cat("\n", paste("--- Starting run for", selVars[1], "x", selVars[2], "- pair", i, "of", nrow(dfBin)), "---\n")

    Data <- select(dat, dfBin[i,])

    colnames(Data)[4] <- "age1"
    colnames(Data)[5] <- "age2"
    Data <- na.omit(Data)

    Data <- as.data.frame(Data, header=T)

    Data[selVars[1]] <- mxFactor(Data[selVars[1]], levels = 0:1, ordered = TRUE)
    Data[selVars[2]] <- mxFactor(Data[selVars[2]], levels = 0:1, ordered = TRUE)

    mzData <- subset(Data, zygos==1)
    dzData <- subset(Data, zygos==2)

    if (nrow(mzData) < min_pairs || nrow(dzData) < min_pairs) {
      cat(sprintf("    WARNING: %s has only %d MZ / %d DZ pairs (< min_pairs = %d). Recording NA row & skipping.\n",
                  selVars[1], nrow(mzData), nrow(dzData), min_pairs))
      uniGN_binary <- rbind(uniGN_binary, data.frame(
        Variable = selVars[1], Label = NA_character_,
        nMZ_pairs = nrow(mzData), nDZ_pairs = nrow(dzData),
        rMZ = NA_real_, rDZ = NA_real_,
        h2 = NA_real_, h2.CI.lower = NA_real_, h2.CI.upper = NA_real_,
        c2 = NA_real_, c2.CI.lower = NA_real_, c2.CI.upper = NA_real_,
        e2 = NA_real_, e2.CI.lower = NA_real_, e2.CI.upper = NA_real_, stringsAsFactors = FALSE))
      next
    }

    nv        <- 1
    nt        <- 2
    nTH      <- 1
    ntv       <- nv*nt

    sVa      <- .3
    sVc      <- .3
    sVe      <- .3
    ThrVals  <- 0

    intercept     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=FALSE, values=0, labels="interC", name="intercept" )
    betaS         <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, labels="betaS", name="bS" )
    betaA         <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values=0, labels="betaA", name="bA" )
    defSex        <- mxMatrix( type="Full", nrow=1, ncol=nt, free=FALSE, labels=c("data.sex1","data.sex2"), name="Sex" )
    defAge        <- mxMatrix( type="Full", nrow=1, ncol=nt, free=FALSE, labels=c("data.age1","data.age2"), name="Age" )

    expMean       <- mxAlgebra( expression = intercept + Sex%x%bS + Age%x%bA , name="expMean" )

    expThr        <- mxMatrix(type="Full", nrow=nTH, ncol=ntv, free = TRUE, values = ThrVals, labels =  paste("th", 1:nTH, sep = ""), name = "expThr")

    covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVa, label="VA11", name="VA" )
    covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVc, label="VC11", name="VC" )
    covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=sVe, label="VE11", name="VE" )
    cons <- mxConstraint(VA+VC+VE ==1, name = "cons")

    covP      <- mxAlgebra( expression= VA+VC+VE, name="V" )
    covMZ     <- mxAlgebra( expression= VA+VC, name="cMZ" )
    covDZ     <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
    expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ),
                                               cbind(t(cMZ), V)), name="expCovMZ" )
    expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ),
                                               cbind(t(cDZ), V)), name="expCovDZ" )

    dataMZ    <- mxData( observed=mzData, type="raw" )
    dataDZ    <- mxData( observed=dzData, type="raw" )

    expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", thresholds = "expThr", dimnames=selVars )
    expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", thresholds = "expThr", dimnames=selVars )
    funML     <- mxFitFunctionML()

    defs      <- list( defAge, defSex )
    pars      <- list( intercept, betaS, betaA, expThr, covA, covC, covE, covP )
    modelMZ   <- mxModel( defs, pars, expMean, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
    modelDZ   <- mxModel( defs, pars, expMean, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
    multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

    rowVC     <- rep('VC',nv)
    colVC     <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
    estVC     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="VarC", dimnames=list(rowVC,colVC) )

    ciACE     <- mxCI( "VarC[1,4:6]" )

    modelACE  <- mxModel("ACEvc", modelMZ, modelDZ, multi, pars, estVC, ciACE, cons )

    fitACE    <- mxRun( modelACE, intervals=T )
    sumACE    <- summary( fitACE )
    mzData[,1] <- as.numeric(mzData[,1])
    mzData[,2] <- as.numeric(mzData[,2])

    dzData[,1] <- as.numeric(dzData[,1])
    dzData[,2] <- as.numeric(dzData[,2])

    nMZ_pairs <- nrow(mzData)
    nDZ_pairs <- nrow(dzData)

    current_nMZ <- nMZ_pairs
    current_nDZ <- nDZ_pairs

    current_rMZ <- fitACE$MZ$cMZ$result[1,1]
    current_rDZ <- fitACE$DZ$cDZ$result[1,1]

    current_h2<-sumACE$CI[1,2]
    current_c2<-sumACE$CI[2,2]
    current_e2<-sumACE$CI[3,2]

    current_h2.CI.lower <- sumACE$CI[1,1]
    current_h2.CI.upper <- sumACE$CI[1,3]

    current_c2.CI.lower <- sumACE$CI[2,1]
    current_c2.CI.upper <- sumACE$CI[2,3]

    current_e2.CI.lower <- sumACE$CI[3,1]
    current_e2.CI.upper <- sumACE$CI[3,3]

    currentVar1_label <- label(dat[selVars[1]])

    uniGN_binary <- rbind(uniGN_binary, data.frame(Variable = selVars[1], Label = currentVar1_label,
                                                   nMZ_pairs = current_nMZ, nDZ_pairs = current_nDZ,
                                                   rMZ = current_rMZ, rDZ = current_rDZ,
                                                   h2=current_h2, h2.CI.lower=current_h2.CI.lower, h2.CI.upper=current_h2.CI.upper,
                                                   c2=current_c2, c2.CI.lower=current_c2.CI.lower, c2.CI.upper=current_c2.CI.upper,
                                                   e2=current_e2, e2.CI.lower=current_e2.CI.lower, e2.CI.upper=current_e2.CI.upper))

  }

  write.csv(uniGN_binary, paste0(outDir, "uniGN_121026_binary.csv"), row.names = FALSE)

  library(writexl)

  uniGN_cont <- read.csv(paste0(outDir, "uniGN_121026_cont.csv"))

  VarPairs_Binary <- list(
    c("ra_level_enrol", "rcqalage1", "A-level enrolment", "Education Achievement & Attainment", "A-level enrolment age 16", "educational attainment"),
    c("rcqhe", "rcqalage1", "university enrolment", "Education Achievement & Attainment", "university enrolment age 18", "educational attainment")
  )

  uniGN_binary$Variable <- sub("1$", "", uniGN_binary$Variable)

  for (p in VarPairs_Binary) {
    idx <- which(uniGN_binary$Variable == p[1])
    if (length(idx) > 0) {
      uniGN_binary$Age_Var_Used[idx]   <- p[2]
      uniGN_binary$Label[idx]          <- p[3]
      uniGN_binary$Category[idx]       <- p[4]
      uniGN_binary$Specific_Label[idx] <- p[5]
      uniGN_binary$Short_Label[idx]    <- p[6]
    }
  }

  missing_cols <- setdiff(names(uniGN_cont), names(uniGN_binary))
  uniGN_binary[missing_cols] <- NA

  uniGN_binary <- uniGN_binary[, names(uniGN_cont)]

  uniGN_all <- rbind(uniGN_cont, uniGN_binary)
  uniGN_all <- unique(uniGN_all)

  uniGN_all$Age_Num <- as.numeric(gsub("\\D*(\\d+).*", "\\1", uniGN_all$Specific_Label))

  cat_levels <- c("Cognitive Abilities", "Education Achievement & Attainment", "Height & BMI", "Psychopathology")
  uniGN_all$Category <- factor(uniGN_all$Category, levels = c(cat_levels, setdiff(unique(uniGN_all$Category), cat_levels)))

  uniGN_all <- uniGN_all[order(uniGN_all$Category, uniGN_all$Age_Num), ]
  uniGN_all$Age_Num <- NULL

  write.csv(uniGN_all, paste0(outDir, "uniGN_121026.csv"), row.names = FALSE)

  uniGN_all_list[[sample_tag]] <- uniGN_all

  cat("Sample [", sample_tag, "] merged. Total variables:", nrow(uniGN_all), "\n")

}

sheet_tables <- list(
  Whole_Sample   = uniGN_all_list[["whole"]],
  ParentDNA = uniGN_all_list[["parentDNA"]]
)
write_xlsx(sheet_tables, paste0(outFileStem, "uniGN_121026.xlsx"))

cat("\nAll done. Two-sheet workbook written to:\n  ",
    paste0(outFileStem, "uniGN_121026.xlsx"), "\n")

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

source("0GN_Theme.R")

for (sample_tag in names(uniGN_all_list)) {

  outDir <- paste0(outFileStem, sample_dirs[[sample_tag]])
  df <- uniGN_all_list[[sample_tag]]

  plot_data <- df %>%
    select(Specific_Label, Category, h2, c2, e2) %>%
    drop_na(h2, c2, e2) %>%
    mutate(Specific_Label = str_replace_all(Specific_Label, "childhood |adolescence |adulthood ", "")) %>%
    mutate(Specific_Label = factor(Specific_Label, levels = unique(Specific_Label)))

  plot_data_long <- plot_data %>%
    pivot_longer(
      cols = c(h2, c2, e2),
      names_to = "Component",
      values_to = "Variance"
    ) %>%
    mutate(
      Component = case_when(
        Component == "h2" ~ "A",
        Component == "c2" ~ "C",
        Component == "e2" ~ "E"
      ),
      Component = factor(Component, levels = c("A", "C", "E"))
    )

  p <- ggplot(plot_data_long, aes(x = Specific_Label, y = Variance, fill = Component)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    facet_wrap(~ Category, scales = "free_x", ncol = 2) +
    scale_fill_manual(
      values = c("A" = "#8A9BAE", "C" = "#1F7A8C", "E" = "#E4E2DD")
    ) +
    labs(
      x = NULL,
      y = "Variance explained"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16, margin = margin(r = 12)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 16),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(t = 8, r = 12, b = 8, l = 8),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.key.size = unit(0.6, "cm")
    )

  ggsave(paste0(outFileStem, "ACE_Variance_Plot_", sample_tag, "_121026.pdf"),
         plot = p, width = 16, height = 11, device = "pdf")

}

construct_map_ace <- c(
  "g age 7"  = "General cognitive ability (g)", "g age 9"  = "General cognitive ability (g)",
  "g age 10" = "General cognitive ability (g)", "g age 12" = "General cognitive ability (g)",
  "g age 16" = "General cognitive ability (g)", "g age 25" = "General cognitive ability (g)",
  "verbal ability age 7"  = "Verbal ability", "verbal ability age 9"  = "Verbal ability",
  "verbal ability age 10" = "Verbal ability", "verbal ability age 12" = "Verbal ability",
  "verbal ability age 16" = "Verbal ability", "verbal ability age 25" = "Verbal ability",
  "nonverbal ability age 7"  = "Nonverbal ability", "nonverbal ability age 9"  = "Nonverbal ability",
  "nonverbal ability age 10" = "Nonverbal ability", "nonverbal ability age 12" = "Nonverbal ability",
  "nonverbal ability age 16" = "Nonverbal ability", "nonverbal ability age 25" = "Nonverbal ability",
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
  "height age 7"  = "Height", "height age 12" = "Height", "height age 14" = "Height",
  "height age 16" = "Height", "height age 22" = "Height", "height age 26" = "Height",
  "BMI age 7"  = "BMI", "BMI age 12" = "BMI", "BMI age 14" = "BMI",
  "BMI age 16" = "BMI", "BMI age 22" = "BMI", "BMI age 26" = "BMI"
)
construct_order_ace8 <- c("Educational attainment", "Educational achievement",
                          "General cognitive ability (g)", "Verbal ability", "Nonverbal ability",
                          "Psychopathology", "Height", "BMI")

for (sample_tag in names(uniGN_all_list)) {

  df <- uniGN_all_list[[sample_tag]]

  ace_trait8 <- df %>%
    select(Specific_Label, Category, h2, c2, e2) %>%
    drop_na(h2, c2, e2) %>%
    mutate(Specific_Label = str_replace_all(Specific_Label, "childhood |adolescence |adulthood ", ""),
           Category       = as.character(Category))

  ace_agg8 <- ace_trait8 %>%
    mutate(Construct = ifelse(Category == "Psychopathology", "Psychopathology",
                              unname(construct_map_ace[Specific_Label])))
  if (any(is.na(ace_agg8$Construct)))
    stop("[", sample_tag, "] traits without a construct: ",
         paste(ace_agg8$Specific_Label[is.na(ace_agg8$Construct)], collapse = ", "))

  ace_agg8 <- ace_agg8 %>%
    group_by(Construct) %>%
    summarise(n_traits = n(), h2 = mean(h2), c2 = mean(c2), e2 = mean(e2), .groups = "drop")
  if (!setequal(ace_agg8$Construct, construct_order_ace8))
    stop("[", sample_tag, "] construct_order_ace8 does not match the constructs in the data")

  ace_agg8_long <- ace_agg8 %>%
    pivot_longer(cols = c(h2, c2, e2), names_to = "Component", values_to = "Variance") %>%
    mutate(Component = factor(recode(Component, h2 = "A", c2 = "C", e2 = "E"),
                              levels = c("A", "C", "E")),
           Construct = factor(Construct, levels = rev(construct_order_ace8)))

  p_agg8 <- ggplot(ace_agg8_long, aes(x = Variance, y = Construct, fill = Component)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.5,
             colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("A" = "#8A9BAE", "C" = "#1F7A8C", "E" = "#E4E2DD"), name = NULL) +
    scale_x_continuous(breaks = BRK_BETA) +
    labs(x = "Variance explained", y = NULL,
         title = "Twin ACE estimates, averaged across ages") +
    theme_gn(grid = "x") +
    guides(fill = guide_legend(nrow = 1))
  ggsave(paste0(outFileStem, "ACE_Variance_Bar_Agg_", sample_tag, "_121026.pdf"),
         plot = p_agg8, width = 16, height = 3.0 + nrow(ace_agg8) * 0.55, device = "pdf")

}
