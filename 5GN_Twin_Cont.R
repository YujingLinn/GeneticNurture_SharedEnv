library(OpenMx)

sourceFileStem <- "GN_Data/"
outFileStem    <- "GN_Results/TwinResults/"

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

source("0GN_VarPairs.R")

binary_traits <- c("ra_level_enrol", "rcqhe")
VarPairs <- Filter(function(v) !(v[1] %in% binary_traits), VarPairs)

for (sample_tag in names(sample_list)) {

  dat    <- sample_list[[sample_tag]]
  outDir <- paste0(outFileStem, sample_dirs[[sample_tag]])
  cat("\n\n########## SAMPLE:", sample_tag, "(N rows =", nrow(dat), ") ##########\n\n")

  uniGN_cont <- data.frame(Variable = character(), Age_Var_Used = character(),
                           Label = character(),
                           Category = character(),
                           Specific_Label = character(),
                           Short_Label = character(),
                           nMZ_pairs = numeric(), nDZ_pairs = numeric(),
                           rMZ = numeric(), rDZ = numeric(),
                           h2=numeric(), h2.CI.lower=numeric(), h2.CI.upper=numeric(),
                           c2=numeric(), c2.CI.lower=numeric(), c2.CI.upper=numeric(),
                           e2=numeric(), e2.CI.lower=numeric(), e2.CI.upper=numeric(),
                           stringsAsFactors = FALSE)

  require(OpenMx)
  mxOption(NULL, "Default optimizer", "SLSQP")

  for (i in 1:length(VarPairs)) {
    current_pair <- VarPairs[[i]]
    currentVar <- current_pair[1]

    age_var_T1 <- current_pair[3]

    if (grepl("1$", age_var_T1)) {
      age_var_T2 <- sub("1$", "2", age_var_T1)
    } else {
      age_var_T2 <- age_var_T1
    }

    meta_label      <- current_pair[4]
    meta_category   <- current_pair[5]
    meta_specific   <- current_pair[6]
    meta_short      <- current_pair[7]

    cat("\n", paste("--- Starting run for", currentVar, "- pair", i, "of", length(VarPairs)), "---\n")
    cat(paste("    Adjusting for Sex and Age: T1(", age_var_T1, "), T2(", age_var_T2, ")\n"))

    var1 <- paste0(currentVar, "1")
    var2 <- paste0(currentVar, "2")

    sex1 <- "sex1"
    sex2 <- "sex2"

    required_cols <- c(var1, var2, sex1, sex2, age_var_T1, age_var_T2, "zygos")
    if (!all(required_cols %in% names(dat))) {
      missing <- required_cols[!required_cols %in% names(dat)]
      cat("    ERROR: Missing columns in dataset:", paste(missing, collapse=", "), "- Recording NA row & skipping.\n")
      uniGN_cont <- rbind(uniGN_cont, data.frame(
        Variable = currentVar, Age_Var_Used = age_var_T1, Label = meta_label,
        Category = meta_category, Specific_Label = meta_specific, Short_Label = meta_short,
        nMZ_pairs = NA_real_, nDZ_pairs = NA_real_, rMZ = NA_real_, rDZ = NA_real_,
        h2 = NA_real_, h2.CI.lower = NA_real_, h2.CI.upper = NA_real_,
        c2 = NA_real_, c2.CI.lower = NA_real_, c2.CI.upper = NA_real_,
        e2 = NA_real_, e2.CI.lower = NA_real_, e2.CI.upper = NA_real_, stringsAsFactors = FALSE))
      next
    }

    work_dat <- dat[, required_cols]

    res_mod_1 <- lm(as.formula(paste(var1, "~", age_var_T1, "+", sex1)),
                    data = work_dat, na.action = na.exclude)
    res_mod_2 <- lm(as.formula(paste(var2, "~", age_var_T2, "+", sex2)),
                    data = work_dat, na.action = na.exclude)

    work_dat[[var1]] <- as.numeric(scale(residuals(res_mod_1)))
    work_dat[[var2]] <- as.numeric(scale(residuals(res_mod_2)))

    nv			<- 1
    ntv			<- nv*2

    selVars <- c(var1, var2)

    mzData		<- subset( work_dat, zygos == 1, selVars )
    dzData		<- subset( work_dat, zygos == 2, selVars )

    mzData <- na.omit(mzData)
    dzData <- na.omit(dzData)

    if (nrow(mzData) < min_pairs || nrow(dzData) < min_pairs) {
      cat(sprintf("    WARNING: only %d MZ / %d DZ pairs (< min_pairs = %d). Recording NA row & skipping.\n",
                  nrow(mzData), nrow(dzData), min_pairs))
      uniGN_cont <- rbind(uniGN_cont, data.frame(
        Variable = currentVar, Age_Var_Used = age_var_T1, Label = meta_label,
        Category = meta_category, Specific_Label = meta_specific, Short_Label = meta_short,
        nMZ_pairs = nrow(mzData), nDZ_pairs = nrow(dzData), rMZ = NA_real_, rDZ = NA_real_,
        h2 = NA_real_, h2.CI.lower = NA_real_, h2.CI.upper = NA_real_,
        c2 = NA_real_, c2.CI.lower = NA_real_, c2.CI.upper = NA_real_,
        e2 = NA_real_, e2.CI.lower = NA_real_, e2.CI.upper = NA_real_, stringsAsFactors = FALSE))
      next
    }

    mzMeanStart	<- mean(mzData[, selVars[1]], na.rm=T)
    dzMeanStart	<- mean(dzData[, selVars[1]], na.rm=T)

    mzCovStart	<- cov(mzData[, selVars], use="complete")[1,2]
    dzCovStart	<- cov(dzData[, selVars], use="complete")[1,2]

    rMZ			<- cor(mzData,use="complete")
    rDZ			<- cor(dzData,use="complete")

    fA			<- 2*(rMZ[2,1] - rDZ[2,1])
    fC			<- rMZ[2,1] - fA
    fE			<- 1 - rMZ[2,1]

    if (fA < 0) { fA = 0 } else if (fA >= 0.95) { fA = 0.9 }
    if (fC < 0) { fC = 0 } else if (fC >= 0.95) { fC = 0.9 }
    if (fE < 0) { fE = 0 } else if (fE >= 0.95) { fE = 0.9 }

    V			<- var(work_dat[[var1]], na.rm=T)

    aStart		<- sqrt(fA * V)
    cStart		<- sqrt(fC * V)
    eStart		<- sqrt(fE * V)

    mzMeanStart = mzMeanStart + (mzMeanStart / 10)
    dzMeanStart = dzMeanStart + (dzMeanStart / 10)
    mzCovStart = mzCovStart + (mzCovStart / 10)
    dzCovStart = dzCovStart + (dzCovStart / 10)
    aStart = aStart + 0.05
    cStart = cStart + 0.05
    eStart = eStart + 0.05

    CholMZ		<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=mzCovStart, name="lowMZ" )
    CholDZ		<-mxMatrix( type="Lower", nrow=ntv, ncol=ntv, free=T, values=dzCovStart, name="lowDZ" )
    MZCov		<-mxAlgebra( expression=lowMZ %*% t(lowMZ), name="expCovMZ" )
    DZCov		<-mxAlgebra( expression=lowDZ %*% t(lowDZ), name="expCovDZ" )

    MZMeans		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=mzMeanStart, labels=c("Mmz1","Mmz2"), name="expMeanMZ" )
    DZMeans		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=dzMeanStart, labels=c("Mdz1","Mdz2"), name="expMeanDZ" )

    matI		<-mxMatrix( type="Iden", nrow=ntv, ncol=ntv, name="I")
    MZcor		<-mxAlgebra( expression= solve(sqrt(I*expCovMZ)) %*% expCovMZ %*% solve(sqrt(I*expCovMZ)), name="rMZ" )
    DZcor		<-mxAlgebra( expression= solve(sqrt(I*expCovDZ)) %*% expCovDZ %*% solve(sqrt(I*expCovDZ)), name="rDZ" )

    dataMZ		<- mxData( observed=mzData, type="raw" )
    dataDZ		<- mxData( observed=dzData, type="raw" )

    objMZ		<- mxExpectationNormal( covariance="expCovMZ", means="expMeanMZ", dimnames=selVars)
    objDZ		<- mxExpectationNormal( covariance="expCovDZ", means="expMeanDZ", dimnames=selVars)

    fitFunction	<- mxFitFunctionML()

    modelMZ		<- mxModel( CholMZ, MZCov, MZMeans, matI, MZcor, dataMZ, objMZ, fitFunction, name="MZ")
    modelDZ		<- mxModel( CholDZ, DZCov, DZMeans, matI, DZcor, dataDZ, objDZ, fitFunction, name="DZ")
    minus2ll	<- mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
    obj			<- mxFitFunctionAlgebra( "m2LL" )
    Conf		<- mxCI (c ('MZ.rMZ[2,1]', 'DZ.rDZ[2,1]') )
    SatModel	<- mxModel( "Sat", modelMZ, modelDZ, minus2ll, obj, Conf)

    SatFit		<- mxRun(SatModel, intervals=T)

    meanG		<- mxMatrix( type="Full", nrow=1, ncol=ntv, free=T, values=mzMeanStart, label="M", name="expMean" )

    pathA		<-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=aStart, label="a11", name="a" )
    pathC		<-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=cStart, label="c11", name="c" )
    pathE		<-mxMatrix( type="Lower", nrow=nv, ncol=nv, free=T, values=eStart, label="e11", name="e" )

    covA		<-mxAlgebra( expression=a %*% t(a), name="A" )
    covC		<-mxAlgebra( expression=c %*% t(c), name="C" )
    covE		<-mxAlgebra( expression=e %*% t(e), name="E" )

    covP		<-mxAlgebra( expression=A+C+E, name="V" )

    Id			<-mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
    SD			<-mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

    stA			<-mxAlgebra( expression=a%*%iSD, name="sta")
    stC			<-mxAlgebra( expression=c%*%iSD, name="stc")
    stE			<-mxAlgebra( expression=e%*%iSD, name="ste")

    h2			<-mxAlgebra( expression=A/V, name="h2")
    c2			<-mxAlgebra( expression=C/V, name="c2")
    e2			<-mxAlgebra( expression=E/V, name="e2")

    estVars <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="Est",
                          dimnames=list(currentVar, c('A','C','E','h2','c2','e2')))

    covMZ		<-mxAlgebra( expression= rbind( cbind(A+C+E,A+C), cbind(A+C,A+C+E)  ), 		name="expCovMZ" )
    covDZ		<-mxAlgebra( expression= rbind( cbind(A+C+E,0.5%x%A+C),cbind(0.5%x%A+C,A+C+E)),	name="expCovDZ" )

    dataMZ		<-mxData( observed=mzData, type="raw" )
    dataDZ		<-mxData( observed=dzData, type="raw" )

    objMZ		<- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars)
    objDZ		<- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars)

    fitFunction	<- mxFitFunctionML()

    params		<-	list(meanG, pathA, pathC, pathE, covA, covC, covE, covP, Id, SD, stA, stC, stE, h2, c2, e2, estVars, fitFunction )
    modelMZ		<-mxModel( params, covMZ, dataMZ, objMZ, name="MZ" )
    modelDZ		<-mxModel( params, covDZ, dataDZ, objDZ, name="DZ" )
    minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
    obj			<-mxFitFunctionAlgebra( "m2LL" )
    Conf		<-mxCI (c ('MZ.h2[1,1]', 'MZ.c2[1,1]', 'MZ.e2[1,1]') )
    ACEModel   	<-mxModel( "ACE", params, modelMZ, modelDZ, minus2ll, obj, Conf)

    ACEFit		<-mxRun(ACEModel, intervals=T)
    (ACESum		<-summary(ACEFit))
    ACE_Path<-round(ACEFit$Est@result,4)

    current_nMZ <- nrow(mzData)
    current_nDZ <- nrow(dzData)

    current_rMZ <- SatFit$MZ$rMZ$result[1,2]
    current_rDZ <- SatFit$DZ$rDZ$result[1,2]

    current_h2<-ACE_Path[1,4]
    current_c2<-ACE_Path[1,5]
    current_e2<-ACE_Path[1,6]

    current_h2.CI.lower <- ACESum$CI[1,1]
    current_h2.CI.upper <- ACESum$CI[1,3]
    current_c2.CI.lower <- ACESum$CI[2,1]
    current_c2.CI.upper <- ACESum$CI[2,3]
    current_e2.CI.lower <- ACESum$CI[3,1]
    current_e2.CI.upper <- ACESum$CI[3,3]

    uniGN_cont <- rbind(uniGN_cont, data.frame(
      Variable = currentVar,
      Age_Var_Used = age_var_T1,
      Label = meta_label,
      Category = meta_category,
      Specific_Label = meta_specific,
      Short_Label = meta_short,
      nMZ_pairs = current_nMZ,
      nDZ_pairs = current_nDZ,
      rMZ = current_rMZ,
      rDZ = current_rDZ,
      h2 = current_h2, h2.CI.lower = current_h2.CI.lower, h2.CI.upper = current_h2.CI.upper,
      c2 = current_c2, c2.CI.lower = current_c2.CI.lower, c2.CI.upper = current_c2.CI.upper,
      e2 = current_e2, e2.CI.lower = current_e2.CI.lower, e2.CI.upper = current_e2.CI.upper
    ))

  }

  write.csv(uniGN_cont, paste0(outDir, "uniGN_121026_cont.csv"), row.names = FALSE)
}

cat("All done (continuous, both samples)! \n")
