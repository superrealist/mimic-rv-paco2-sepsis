# Shared paths and fixed study settings. Run scripts from the repository root.
CODE <- normalizePath(Sys.getenv('MIMIC_CODE_ROOT','.'),winslash='/',mustWork=TRUE)
DATA <- normalizePath(Sys.getenv('MIMIC_DATA_DIR',file.path(CODE,'data')),winslash='/',mustWork=FALSE)
OUT <- normalizePath(Sys.getenv('MIMIC_OUTPUT_DIR',file.path(CODE,'outputs')),winslash='/',mustWork=FALSE)
INPUT <- Sys.getenv('MIMIC_ANALYSIS_INPUT',file.path(DATA,'analysis_input.csv'))
WORKERS <- as.integer(Sys.getenv('MIMIC_WORKERS','4'))
stopifnot(is.finite(WORKERS),WORKERS>=1L)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
FUN <- function(name) file.path(CODE,'R/functions',name)
KN <- list(age=c(42,66,84),bmi=c(21.484776785714285,29.038986884611216,41.5497448979592),
           map=c(60,74,93),pf=c(103.75,210,375))
COVARIATES <- 'age+agenl+bmi+bminl+map+mapnl+pf+pfnl+llac+lmv+sofa+peep+lvent+sex+copd+chf+ph+cad+ckd'
FORM_PRIMARY <- paste0('y~rv*(x+xnl)+',COVARIATES)
FORM_CATEGORY <- paste0('y~rv*(co2lo+co2hi)+',COVARIATES)
FORM_CATEGORY_IMPUTATION <- paste0(FORM_CATEGORY,'+hist')
M_PRIMARY <- 100L
M_CATEGORY <- 100L
M_SUPPORT <- 40L
ITERATIONS <- 40L
SUPPORT_ENDPOINTS <- c('vpfd7','vfd7','crrtfd7','vpfd28','vfd28','crrtfd28')
make_dirs <- function(root,names) for(n in names) dir.create(file.path(root,n),recursive=TRUE,showWarnings=FALSE)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
