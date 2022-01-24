#!/usr/bin/env Rscript

##### NOTE THAT 3 VARIABLES HAVE NO DATA, WE'RE ONLY LOOKING AT 21 MAX

library(tidyverse)
library(ComplexUpset)

logged = readRDS('logged_cohort1_12_to_24_var_voe_T2D.rds')
logged_processed = map(logged, function(x) x[[2]] %>% mutate(significant_initial = if_else(BY<=0.05,1,0)) %>% select(feature,significant_initial)) %>% reduce(inner_join,by='feature')
colnames(logged_processed) = c('feature','var_12','var_15','var_18','var_21','var_24')
logged_processed = logged_processed %>% select(-var_12,-feature)

pdf(paste('12up_overall_sigonce.pdf',sep=''),width=12,height=9)
upset(logged_processed %>% as.data.frame,colnames(logged_processed),keep_empty_groups=TRUE)
dev.off()

### pval sig at least once plots
logged_processed = map(logged, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_processed) = c('feature','var_12','var_15','var_18','var_21','var_24')
logged_processed = logged_processed %>% select(-var_12,-feature)

pdf(paste('12up_overall_pvalsigonce.pdf',sep=''),width=12,height=9)
upset(logged_processed %>% as.data.frame,colnames(logged_processed),keep_empty_groups=TRUE)
dev.off()

### FDR sig at least once plot
bonf_equivalent = 0.05/535
logged_processed = map(logged, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_processed) = c('feature','var_12','var_15','var_18','var_21','var_24')
logged_processed = logged_processed %>% select(-var_12,-feature)

pdf(paste('12up_fdrsig.pdf',sep=''),width=12,height=9)
upset(logged_processed %>% as.data.frame,colnames(logged_processed),keep_empty_groups=TRUE, sort_sets=FALSE)
dev.off()
