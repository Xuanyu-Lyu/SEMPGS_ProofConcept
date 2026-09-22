library('polycor') # For polyserial(). Can comment this out if not using the r2_polyserial() function.
suppressPackageStartupMessages(library('OpenMx'))

# This file contains lots of functions for calculating r2 (for variance a PGS or other continuous variable explains in a phenotype)
# For continuous phenos:
#   * r2_linreg() will get the r2 for a linear regression model. For PGS and continuous phenotype WITHOUT
#     covars, or for use by other functions when calculating incremental r2s.
#   * r2_linreg_incremental() will get the incremental r2 from linear regression using the standard formula for incremental r2 (r2_full - r2_reduced)
#   * r2_linreg_partial() will calculate (r2_full - r2_reduced)/(1 - r2_reduced), producing an r2 similar to a squared partial correlation.
#     The full model is predictor + covars, the reduced model is covars-only.
#     THIS IS WHAT MATT WANTED FOR CONTINUOUS PHENOS WHEN USING COVARIATES.
# Main functions for binary phenos:
#   * r2_polyserial() will get a quick *univariate* r2 for a binary phenotype without needing OpenMx
#   * r2_omx_univar() will get a *univariate* r2 for a binary phenotype using OpenMx
#   * r2_omx_partial() will get an r2 for a binary phenotype using OpenMx and accounting for covariates
# Additional functions for binary phenos using probit regression (experimental, recommended to stick to the other functions above instead):
#   * r2_probit() get total r2 explained by predictor variables. For univariate use or for use by partial/incremental functions.
#   * r2_probit_incremental() will get the incremental r2 from probit regression using the standard formula for incremental r2 (r2_full - r2_reduced)
#   * r2_probit_partial() will calculate (r2_full - r2_reduced)/(1 - r2_reduced), producing an r2 similar to a squared partial correlation.

# Guidelines for use for PGS r2:
# * If you don't have covariates, use r2_linreg() and r2_omx_univar()
# * If you have covariates, use r2_linreg_partial() and r2_omx_partial()


########## Functions for continuous phenotypes ##########

# Get r2 for a linear regression model
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predvars: Column name for predictor variable(s) (can be just the PGS, or a vector of PGS + covars)
r2_linreg <- function(samp,pheno,predvars) {
  reg_formula <- as.formula(paste0(pheno,' ~ ',paste(predvars,collapse=' + ')))
  r2 <- summary(lm(reg_formula,data=samp))$r.squared
  return(r2)
}

# Incremental r2 using linear regression. r2_full - r2_covars
# Full model is predictor + covars, reduced model is covars-only.
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predictor: Column name for predictor (eg. the PGS)
#    covars: Vector of column names for the covariates
r2_linreg_incremental <- function(samp,pheno,predictor,covars) {
  r2_reduced <- r2_linreg(samp,pheno,covars)
  r2_full <- r2_linreg(samp,pheno,c(predictor,covars))
  r2 <- r2_full - r2_reduced
  return(r2)
}

# Partial r2 using linear regression. Denominator is (1 - r2_reduced) to mimic squared partial correlation.
# Full model is predictor + covars, reduced model is covars-only.
# THIS IS WHAT MATT WANTED FOR CONTINUOUS PHENOS WHEN USING COVARIATES.
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predictor: Column name for predictor (eg. the PGS)
#    covars: Vector of column names for the covariates
r2_linreg_partial <- function(samp,pheno,predictor,covars) {
  r2_reduced <- r2_linreg(samp,pheno,covars)
  r2_full <- r2_linreg(samp,pheno,c(predictor,covars))
  r2 <- (r2_full - r2_reduced)/(1 - r2_reduced)
  return(r2)
}


########## Main functions for binary phenotypes ##########

# Get r2 for a binary phenotype based on squaring a polyserial correlation. *NO* covariates or partial correlation.
# Polyserial treats the binary or categorical variable as ordinal and looks at correlation with underlying liability
# Requires the polyserial() function from the 'polycor' library.
#    samp: A data.frame containing the phenotype and predictor of interest (PGS)
#    pheno: The name of the column containing the binary phenotype
#    predictor: Column name for predictor (eg. the PGS)
#    ML: (T/F) Whether to use maximum likelihood when running polyserial(). More accurate but slower.
r2_polyserial <- function(samp,pheno,predictor,ML=T) {
  r2 <- polyserial(samp[,predictor],samp[,pheno],ML=ML)^2
  return(r2)
}

# Univariate r2 with OpenMx
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predictor: Column name for predictor (eg. the PGS)
#    scale_predictor: T/F, make predictor have a mean of 0 and variance of 1. Not absolutely necessary because we use variance to convert covariance to correlation, but helps with model fitting/start values.
r2_omx_univar <- function(samp,pheno,predictor,scale_predictor=T) {
  indat <- samp[,c(pheno,predictor)]
  if(is.logical(indat[,pheno])) { # If currently T/F instead of 0/1, convert to 0/1
    indat[,pheno] <- 1*indat[,pheno]
  }
  stopifnot(all(indat[,pheno] %in% c(0,1,NA))) # Function not meant for contin phenos
  indat[,pheno] <- mxFactor(as.vector(indat[,pheno]),levels=c(0,1))
  if(scale_predictor) { # Make predictor have a mean of 0 and variance of 1
    indat[,predictor] <- c(scale(indat[,predictor],center=T,scale=T))
  }
  # Input matrices for the model
  covmat <- mxMatrix('Symm',2,2,byrow=T,name='covmat',
                      values=c(1,0.2,
                               0.2,1),
                      free=c(T,T,
                             T,F),
                      labels=c('var_contin','offdiag',
                               'offdiag','fixed_latent_variance'),
                      dimnames=list(c(predictor,pheno),c(predictor,pheno)))
  meanmat <- mxMatrix('Full',nrow=1,ncol=2,values=c(0,0),free=c(T,F),dimnames=list(NULL,c(predictor,pheno)),name='meanmat')
  z_thresh <- qnorm(mean(indat[,pheno]=='0',na.rm=T)) # Estimate z-score corresponding to liability threshold, based on prevalence
  threshmat <- mxMatrix('Full',1,1,free=T,values=z_thresh,dimnames=list(pheno,NULL),name='threshmat') # 1x1 matrix. Matt's comments said no column names to avoid name mismatches.
  # Define how model expectations are calculated. Matt's comments for threshnames= said "explicitly state which var has thresholds"
  expectation <- mxExpectationNormal(covariance='covmat',means='meanmat',dimnames=c(predictor,pheno),threshold='threshmat',threshnames=pheno)
  # Put the pieces together into a model, then fit/run the model
  model <- mxModel('polyserial_ML',mxData(indat, type = "raw"),covmat,meanmat,threshmat,expectation,mxFitFunctionML())
  fit <- mxRun(model)
  # Get results
  covariance <- mxEval(offdiag, fit)
  # covariance_se  <- mxSE(covariance[1,1], fit) # Not currently used, but worth keeping example
  var_predictor <- mxEval(var_contin, fit)
  rho_hat <- covariance/sqrt(var_predictor) # cor(x,y) = cov(x,y)/(se(x)*se(y)). We force variance of latent pheno liability to 1, so only need variance/SE of PGS here. Must be variance AFTER removing covars.
  r2 <- rho_hat^2
  return(r2)
}

# Partial r2 with OpenMx and covars
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predictor: Column name for predictor (eg. the PGS)
#    covars: Vector of column names for covariates
#    scale_predictor: T/F, make predictor have a mean of 0 and variance of 1. Not absolutely necessary because we use variance to convert covariance to correlation, but helps with model fitting/start values.
r2_omx_partial <- function(samp,pheno,predictor,covars,scale_predictor=T) {
  indat <- samp[,c(pheno,predictor,covars)]
  if(is.logical(indat[,pheno])) { # If currently T/F instead of 0/1, convert to 0/1
    indat[,pheno] <- 1*indat[,pheno]
  }
  stopifnot(all(indat[,pheno] %in% c(0,1,NA))) # Function not meant for contin phenos
  indat[,pheno] <- mxFactor(as.vector(indat[,pheno]),levels=c(0,1))
  if(scale_predictor) { # Make predictor have a mean of 0 and variance of 1.
    indat[,predictor] <- c(scale(indat[,predictor],center=T,scale=T))
  }
  # Matt's comment: Scalar matrices (NO labels to avoid name collisions) 
  mu    <- mxMatrix("Full", 1, 1, free=TRUE,  values=0,   name="mu") # Mean of continuous predictor (the PGS) after regressing out covariates
  betas <- lapply(seq_along(covars), function(i) { mxMatrix(type="Full", nrow=1, ncol=1, values=0, free=TRUE, name=paste0("b", i)) }) # Betas for covariates
  t0    <- mxMatrix("Full", 1, 1, free=TRUE,  values=0,   name="t0") # t0 for the threshold (threshold after regressing out covariates)
  gs <- lapply(seq_along(covars), function(i) { mxMatrix(type="Full", nrow=1, ncol=1, values=0, free=TRUE, name=paste0("g", i)) }) # Coefficients for covariates in the threshold mxAlgebra
  vpredictor    <- mxMatrix("Full", 1, 1, free=TRUE,  values=1,   name="vpredictor") # Variance for the continuous predictor (the PGS)
  offdiag   <- mxMatrix("Full", 1, 1, free=TRUE,  values=.2,  name="offdiag") # Covariance between the continuous predictor of interest and the underlying liability for the phenotype
  Unit1 <- mxMatrix("Full", 1, 1, free=FALSE, values=1,   name="Unit1") # Variance for the underlying liability for the phenotype (fixed to 1)
  # Means
  mean_formula <- paste0('mu + ',paste(paste0('b',1:length(covars),'*data.',covars),collapse=' + ')) # mu + b1*data.sex + b2*data.covar1, etc.
  MEANcont <- mxAlgebraFromString(mean_formula, name="MEANcont")
  meanalg <- mxAlgebra(cbind(MEANcont, 0),dimnames=list(NULL, c(predictor,pheno)),name="meanalg") # This links MEANcont to the observed mean for the predictor, and sets the mean for the underlying liability for the pheno to 0
  # Covariance matrix for continuous predictor and liability for binary pheno
  covalg <- mxAlgebra(rbind(cbind(vpredictor, offdiag),
                                  cbind(offdiag, Unit1)),
                            dimnames=list(c(predictor,pheno),c(predictor,pheno)),name='covalg')
  # Threshold as regression: tau = t0 + g1*cov1 + g2*cov2
  thresh_formula <- paste0('t0 + ',paste(paste0('g',1:length(covars),'*data.',covars),collapse=' + ')) # t0 + g1*data.sex + g2*data.covar1, etc.
  threshalg <- mxAlgebraFromString(thresh_formula,dimnames=list(pheno, NULL),name='threshalg')
  # Define how model expectations are calculated. Matt's comments for threshnames= said "explicitly state which var has thresholds"
  expectation <- mxExpectationNormal(covariance='covalg',means='meanalg',dimnames=c(predictor,pheno),thresholds='threshalg',threshnames=pheno, jointConditionOn = "continuous")
  # Put the pieces together into a model, then fit/run the model
  model <- mxModel("partial_polyserial_defvars",mxData(indat, type = "raw"),
                   mu,betas,t0,gs,vpredictor,offdiag,Unit1,MEANcont,meanalg,covalg,threshalg,expectation,mxFitFunctionML())
  fit <- mxRun(model)
  # Get results
  covariance <- mxEval(offdiag[1,1], fit)
  # covariance_se  <- mxSE(covariance[1,1], fit) # Not currently used, but worth keeping example
  var_predictor <- mxEval(vpredictor[1,1], fit)
  rho_hat <- covariance/sqrt(var_predictor) # cor(x,y) = cov(x,y)/(se(x)*se(y)). We force variance of latent pheno liability to 1, so only need variance/SE of PGS here. Must be variance AFTER removing covars.
  r2 <- rho_hat^2
  return(r2)
}


########## Functions that use probit regression to do r2 for binary phenotypes ##########
# These functions are an optional alternative to the omx functions, but are still experimental
# In initial testing, univariate r2_probit() beat r2_omx_univar() for which matched more closely to r2_polyserial(), but they were both usually pretty close
# BUT, r2_probit_partial() can produce negative values when r2 is close to 0, so I think the OpenMx version is probably preferable in most cases

# Use probit regression and McKelvey-Zavoina R-squared to get total variance predvars explain in latent liability to a binary pheno
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predvars: Column name for predictor variable(s) (can be just the PGS, or a vector of PGS + covars)
r2_probit <- function(samp,pheno,predvars) {
  samp <- samp[rowSums(is.na(samp[,c(pheno,predvars)]))==0,] # Drop missings so they don't mess with results from predict() below
  reg_formula <- as.formula(paste0(pheno,' ~ ',paste(predvars,collapse=' + ')))
  model <- glm(reg_formula, data=samp, family=binomial(link="probit"))
  # Get predicted liability for each person
  pred_liab <- predict(model,type='link') - coef(model)[1] # for each person, sum of predictors*betas, minus intercepts. type='link' keeps it on the scale of the predictors rather than making it a predicted probability
  var_pred <- var(pred_liab)
  r2 <- var_pred / (var_pred + 1) # Residual variance for probit regression is always 1
  return(r2)
}

# Incremental r2 using probit regression. r2_full - r2_covars
# Full model is predictor + covars, reduced model is covars-only.
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predictor: Column name for predictor (eg. the PGS)
#    covars: Vector of column names for the covariates
r2_probit_incremental <- function(samp,pheno,predictor,covars) {
  r2_reduced <- r2_probit(samp,pheno,covars)
  r2_full <- r2_probit(samp,pheno,c(predictor,covars))
  r2 <- r2_full - r2_reduced
  return(r2)
}

# Partial r2 using probit regression. Denominator is (1 - r2_reduced) to mimic squared partial correlation.
# Full model is predictor + covars, reduced model is covars-only.
#    samp: A data.frame containing pheno and predvars
#    pheno: Column name for the phenotype
#    predictor: Column name for predictor (eg. the PGS)
#    covars: Vector of column names for the covariates
r2_probit_partial <- function(samp,pheno,predictor,covars) {
  r2_reduced <- r2_probit(samp,pheno,covars)
  r2_full <- r2_probit(samp,pheno,c(predictor,covars))
  r2 <- (r2_full - r2_reduced)/(1 - r2_reduced)
  return(r2)
}