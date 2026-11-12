source(file.path(Sys.getenv('MIMIC_CODE_ROOT','.'),'R/config.R'))
O<-file.path(OUT,'contrasts')
source(FUN('pooling.R'));source(FUN('risk_contrasts.R'))
f<-file.path(OUT,'primary/fits/model.rds')
z<-readRDS(f);datasets<-z$data;e<-pi_extract(z$fits)
src<-read.csv(INPUT,check.names=FALSE)
stopifnot(length(datasets)==100L,nrow(datasets[[1]])==2051L,
 identical(z$sourceIDs,src$subject_id),all(src$death28==datasets[[1]]$y))
form<-paste(deparse(formula(z$fits[[1]])),collapse=' ')
calculate_contrasts('MIMIC-IV',datasets,e$Q,e$U,form,'mimic_raw','mimic',e$dfcom,O,
 source_metadata=list(input_sha256=digest::digest(file=INPUT,algo='sha256'),
 model_sha256=digest::digest(file=f,algo='sha256')),workers=WORKERS)
tests<-read.csv(file.path(OUT,'primary/tables/tests.csv'))
write.csv(data.frame(cohort='MIMIC-IV',test=tests$test,df=tests$df1,p_value=tests$p.value,
 statistic=tests$statistic,df_den=tests$df2),file.path(O,'tables/mimic_tests.csv'),row.names=FALSE)
