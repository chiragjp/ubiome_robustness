#!/usr/bin/env Rscript

library(tidyverse)
library(ComplexUpset)
library(pheatmap)

remove_meta_analyzed <- function(data){
	count = 0 
	for(val in data){
		count = count + 1
		if('meta_analyis_output' %in% names(val)){
			data[[count]] = val[names(val)!='initial_association_output']
		}
	}
	return(data)
}

args = commandArgs(trailingOnly=TRUE)

files = read.csv(args[[1]],header=F) %>% unlist %>% unname

clr_3var = readRDS(files[[1]])

ma=FALSE
matruecheck = map(clr_3var, function(x) length(x)) %>% unique
if(6 %in% matruecheck){
	ma=TRUE
}

clr_3var = remove_meta_analyzed(clr_3var)
clr_6var = readRDS(files[[2]])
clr_6var = remove_meta_analyzed(clr_6var)
clr_9var = readRDS(files[[3]])
clr_9var = remove_meta_analyzed(clr_9var)

logged_3var = readRDS(files[[4]])
logged_3var = remove_meta_analyzed(logged_3var)
logged_6var = readRDS(files[[5]])
logged_6var = remove_meta_analyzed(logged_6var)
logged_9var = readRDS(files[[6]])
logged_9var = remove_meta_analyzed(logged_9var)

nonlogged_3var = readRDS(files[[7]])
nonlogged_3var = remove_meta_analyzed(nonlogged_3var)
nonlogged_6var = readRDS(files[[8]])
nonlogged_6var = remove_meta_analyzed(nonlogged_6var)
nonlogged_9var = readRDS(files[[9]])
nonlogged_9var = remove_meta_analyzed(nonlogged_9var)

clr_3var_processed = map(clr_3var, function(x) x[[2]] %>% mutate(significant_initial = if_else(BY<=0.05,1,0)) %>% select(feature,significant_initial)) %>% reduce(inner_join,by='feature')

dset_num = ncol(clr_3var_processed) - 1
clr_cols = c('feature',paste('clr_dataset_',seq(1,dset_num),sep=''))
logged_cols = c('feature',paste('logged_dataset_',seq(1,dset_num),sep=''))
nl_cols = c('feature',paste('nonlogged_dataset_',seq(1,dset_num),sep=''))

if(ma==TRUE){
	dset_num = ncol(clr_3var_processed) - 2
	clr_cols = c('feature',paste('clr_dataset_',seq(1,dset_num),sep=''),'clr_meta_analysis')
	logged_cols = c('feature',paste('logged_dataset_',seq(1,dset_num),sep=''),'logged_meta_analysis')
	nl_cols = c('feature',paste('nonlogged_dataset_',seq(1,dset_num),sep=''),'nonlogged_meta_analysis')
}

colnames(clr_3var_processed) = clr_cols

logged_3var_processed = map(logged_3var, function(x) x[[2]] %>% mutate(significant_initial = if_else(BY<=0.05,1,0)) %>% select(feature,significant_initial)) %>% reduce(inner_join,by='feature')
colnames(logged_3var_processed) = logged_cols

nonlogged_3var_processed = map(nonlogged_3var, function(x) x[[2]] %>% mutate(significant_initial = if_else(BY<=0.05,1,0)) %>% select(feature,significant_initial)) %>% reduce(inner_join,by='feature')
colnames(nonlogged_3var_processed) = nl_cols

var3_all = full_join(clr_3var_processed,logged_3var_processed,by='feature')
var3_all = full_join(var3_all,nonlogged_3var_processed,by='feature') %>% select(-feature)
var3_all[is.na(var3_all)] = 0

pdf(paste(args[[1]],'overall_pvalsig.pdf',sep=''),width=12,height=9)
upset(var3_all %>% as.data.frame,colnames(var3_all),keep_empty_groups=TRUE)
dev.off()

### pval sig at least once plots

dset_num = ncol(clr_3var_processed) - 1
clr_cols = c('feature',paste('clr_3var_dataset_',seq(1,dset_num),sep=''))
logged_cols = c('feature',paste('logged_3var_dataset_',seq(1,dset_num),sep=''))
nl_cols = c('feature',paste('nonlogged_3var_dataset_',seq(1,dset_num),sep=''))

if(ma==TRUE){
	dset_num = ncol(clr_3var_processed) - 2
	clr_cols = c('feature',paste('clr_3var_dataset_',seq(1,dset_num),sep=''),'clr_meta_analysis')
	logged_cols = c('feature',paste('logged_3var_dataset_',seq(1,dset_num),sep=''),'logged_meta_analysis')
	nl_cols = c('feature',paste('nonlogged_3var_dataset_',seq(1,dset_num),sep=''),'nonlogged_meta_analysis')
}

# 3 var
clr_3var_processed = map(clr_3var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(clr_3var_processed) = clr_cols

logged_3var_processed = map(logged_3var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_3var_processed) = logged_cols

nonlogged_3var_processed = map(nonlogged_3var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_3var_processed) = nl_cols

var3_all = full_join(clr_3var_processed,logged_3var_processed,by='feature')
var3_all = full_join(var3_all,nonlogged_3var_processed,by='feature') %>% select(-feature)
var3_all[is.na(var3_all)] = 0

pdf(paste(args[[1]],'3var_pvalsig.pdf',sep=''),width=12,height=9)
upset(var3_all %>% as.data.frame,colnames(var3_all),keep_empty_groups=TRUE, sort_sets=FALSE)
dev.off()

# 6 var

dset_num = ncol(clr_3var_processed) - 1
clr_cols = c('feature',paste('clr_6var_dataset_',seq(1,dset_num),sep=''))
logged_cols = c('feature',paste('logged_6var_dataset_',seq(1,dset_num),sep=''))
nl_cols = c('feature',paste('nonlogged_6var_dataset_',seq(1,dset_num),sep=''))

if(ma==TRUE){
	dset_num = ncol(clr_3var_processed) - 2
	clr_cols = c('feature',paste('clr_6var_dataset_',seq(1,dset_num),sep=''),'clr_meta_analysis')
	logged_cols = c('feature',paste('logged_6var_dataset_',seq(1,dset_num),sep=''),'logged_meta_analysis')
	nl_cols = c('feature',paste('nonlogged_6var_dataset_',seq(1,dset_num),sep=''),'nonlogged_meta_analysis')
}

clr_6var_processed = map(clr_6var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(clr_6var_processed) = clr_cols

logged_6var_processed = map(logged_6var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_6var_processed) = logged_cols

nonlogged_6var_processed = map(nonlogged_6var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_6var_processed) = nl_cols

var6_all = full_join(clr_6var_processed,logged_6var_processed,by='feature')
var6_all = full_join(var6_all,nonlogged_6var_processed,by='feature') %>% select(-feature)
var6_all[is.na(var6_all)] = 0

pdf(paste(args[[1]],'6var_pvalsig.pdf',sep=''),width=12,height=9)
upset(var6_all %>% as.data.frame,colnames(var6_all),keep_empty_groups=TRUE, sort_sets=FALSE)
dev.off()

# 9 var

dset_num = ncol(clr_3var_processed) - 1
clr_cols = c('feature',paste('clr_9var_dataset_',seq(1,dset_num),sep=''))
logged_cols = c('feature',paste('logged_9var_dataset_',seq(1,dset_num),sep=''))
nl_cols = c('feature',paste('nonlogged_9var_dataset_',seq(1,dset_num),sep=''))

if(ma==TRUE){
	dset_num = ncol(clr_3var_processed) - 2
	clr_cols = c('feature',paste('clr_9var_dataset_',seq(1,dset_num),sep=''),'clr_meta_analysis')
	logged_cols = c('feature',paste('logged_9var_dataset_',seq(1,dset_num),sep=''),'logged_meta_analysis')
	nl_cols = c('feature',paste('nonlogged_9var_dataset_',seq(1,dset_num),sep=''),'nonlogged_meta_analysis')
}


clr_9var_processed = map(clr_9var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(clr_9var_processed) = logged_cols

logged_9var_processed = map(logged_9var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_9var_processed) = clr_cols

nonlogged_9var_processed = map(nonlogged_9var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<0.05,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_9var_processed) = nl_cols

var9_all = full_join(clr_9var_processed,logged_9var_processed,by='feature')
var9_all = full_join(var9_all,nonlogged_9var_processed,by='feature') %>% select(-feature)
var9_all[is.na(var9_all)] = 0

pdf(paste(args[[1]],'9var_pvalsig.pdf',sep=''),width=12,height=9)
upset(var9_all %>% as.data.frame,colnames(var9_all),keep_empty_groups=TRUE, sort_sets=FALSE)
dev.off()
### fdr sig at least once plots

# 3 var
dset_num = ncol(clr_3var_processed) - 1
clr_cols = c('feature',paste('clr_3var_dataset_',seq(1,dset_num),sep=''))
logged_cols = c('feature',paste('logged_3var_dataset_',seq(1,dset_num),sep=''))
nl_cols = c('feature',paste('nonlogged_3var_dataset_',seq(1,dset_num),sep=''))

if(ma==TRUE){
	dset_num = ncol(clr_3var_processed) - 2
	clr_cols = c('feature',paste('clr_3var_dataset_',seq(1,dset_num),sep=''),'clr_meta_analysis')
	logged_cols = c('feature',paste('logged_3var_dataset_',seq(1,dset_num),sep=''),'logged_meta_analysis')
	nl_cols = c('feature',paste('nonlogged_3var_dataset_',seq(1,dset_num),sep=''),'nonlogged_meta_analysis')
}

bonf_equivalent = 0.05/535
clr_3var_processed = map(clr_3var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(clr_3var_processed) = clr_cols

logged_3var_processed = map(logged_3var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_3var_processed) = logged_cols

nonlogged_3var_processed = map(nonlogged_3var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_3var_processed) = nl_cols

var3_all = full_join(clr_3var_processed,logged_3var_processed,by='feature')
var3_all = full_join(var3_all,nonlogged_3var_processed,by='feature') %>% select(-feature)
var3_all[is.na(var3_all)] = 0

pdf(paste(args[[1]],'3var_fdrsig.pdf',sep=''),width=12,height=9)
upset(var3_all %>% as.data.frame,colnames(var3_all),keep_empty_groups=TRUE, sort_sets=FALSE)
dev.off()

# 6 var

dset_num = ncol(clr_3var_processed) - 1
clr_cols = c('feature',paste('clr_6var_dataset_',seq(1,dset_num),sep=''))
logged_cols = c('feature',paste('logged_6var_dataset_',seq(1,dset_num),sep=''))
nl_cols = c('feature',paste('nonlogged_6var_dataset_',seq(1,dset_num),sep=''))

if(ma==TRUE){
	dset_num = ncol(clr_3var_processed) - 2
	clr_cols = c('feature',paste('clr_6var_dataset_',seq(1,dset_num),sep=''),'clr_meta_analysis')
	logged_cols = c('feature',paste('logged_6var_dataset_',seq(1,dset_num),sep=''),'logged_meta_analysis')
	nl_cols = c('feature',paste('nonlogged_6var_dataset_',seq(1,dset_num),sep=''),'nonlogged_meta_analysis')
}


clr_6var_processed = map(clr_6var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(clr_6var_processed) = clr_cols

logged_6var_processed = map(logged_6var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_6var_processed) = logged_cols

nonlogged_6var_processed = map(nonlogged_6var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_6var_processed) = nl_cols

var6_all = full_join(clr_6var_processed,logged_6var_processed,by='feature')
var6_all = full_join(var6_all,nonlogged_6var_processed,by='feature') %>% select(-feature)
var6_all[is.na(var6_all)] = 0

pdf(paste(args[[1]],'6var_fdrsig.pdf',sep=''),width=12,height=9)
upset(var6_all %>% as.data.frame,colnames(var6_all),keep_empty_groups=TRUE, sort_sets=FALSE)
dev.off()

# 9 var
dset_num = ncol(clr_3var_processed) - 1
clr_cols = c('feature',paste('clr_9var_dataset_',seq(1,dset_num),sep=''))
logged_cols = c('feature',paste('logged_9var_dataset_',seq(1,dset_num),sep=''))
nl_cols = c('feature',paste('nonlogged_9var_dataset_',seq(1,dset_num),sep=''))

if(ma==TRUE){
	dset_num = ncol(clr_3var_processed) - 2
	clr_cols = c('feature',paste('clr_9var_dataset_',seq(1,dset_num),sep=''),'clr_meta_analysis')
	logged_cols = c('feature',paste('logged_9var_dataset_',seq(1,dset_num),sep=''),'logged_meta_analysis')
	nl_cols = c('feature',paste('nonlogged_9var_dataset_',seq(1,dset_num),sep=''),'nonlogged_meta_analysis')
}


clr_9var_processed = map(clr_9var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(clr_9var_processed) = clr_cols

logged_9var_processed = map(logged_9var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(logged_9var_processed) = logged_cols

nonlogged_9var_processed = map(nonlogged_9var, function(x) x[[5]][[3]]  %>% mutate(pval_sig = if_else(p.value<bonf_equivalent,1,0)) %>% select(dependent_feature,pval_sig) %>%  group_by(dependent_feature) %>% summarize(pval_sig_once = sum(pval_sig)) %>% mutate(pval_sig_once = if_else(pval_sig_once == 0, 0, 1))) %>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_9var_processed) = nl_cols

var9_all = full_join(clr_9var_processed,logged_9var_processed,by='feature')
var9_all = full_join(var9_all,nonlogged_9var_processed,by='feature') %>% select(-feature)
var9_all[is.na(var9_all)] = 0

pdf(paste(args[[1]],'9var_fdrsig.pdf',sep=''),width=12,height=9)
upset(var9_all %>% as.data.frame,colnames(var9_all),keep_empty_groups=TRUE, sort_sets=FALSE)
dev.off()



### heatmap process
# 3 var
clr_3var_processed = map(clr_3var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(clr_3var_processed) = c('feature','clr_3var_dataset_1','clr_3var_dataset_2','clr_3var_dataset_3')

logged_3var_processed = map(logged_3var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(logged_3var_processed) = c('feature','logged_3var_dataset_1','logged_3var_dataset_2','logged_3var_dataset_3')

nonlogged_3var_processed = map(nonlogged_3var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_3var_processed) = c('feature','nonlogged_3var_dataset_1','nonlogged_3var_dataset_2','nonlogged_3var_dataset_3')

var3_all = full_join(clr_3var_processed,logged_3var_processed,by='feature')
var3_all = full_join(var3_all,nonlogged_3var_processed,by='feature')
var3_all[is.na(var3_all)] = 0

# 6 var 
clr_6var_processed = map(clr_6var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(clr_6var_processed) = c('feature','clr_6var_dataset_1','clr_6var_dataset_2','clr_6var_dataset_3')

logged_6var_processed = map(logged_6var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(logged_6var_processed) = c('feature','logged_6var_dataset_1','logged_6var_dataset_2','logged_6var_dataset_3')

nonlogged_6var_processed = map(nonlogged_6var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_6var_processed) = c('feature','nonlogged_6var_dataset_1','nonlogged_6var_dataset_2','nonlogged_6var_dataset_3')

var6_all = full_join(clr_6var_processed,logged_6var_processed,by='feature')
var6_all = full_join(var6_all,nonlogged_6var_processed,by='feature')
var6_all[is.na(var6_all)] = 0

# 9 var
clr_9var_processed = map(clr_9var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(clr_9var_processed) = c('feature','clr_9var_dataset_1','clr_9var_dataset_2','clr_9var_dataset_3')

logged_9var_processed = map(logged_9var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(logged_9var_processed) = c('feature','logged_9var_dataset_1','logged_9var_dataset_2','logged_9var_dataset_3')

nonlogged_9var_processed = map(nonlogged_9var[1:dset_num], function(x) x[[5]][[1]] %>% select(dependent_feature,janus_effect))%>% reduce(inner_join,by='dependent_feature')
colnames(nonlogged_9var_processed) = c('feature','nonlogged_9var_dataset_1','nonlogged_9var_dataset_2','nonlogged_9var_dataset_3')

var9_all = full_join(clr_9var_processed,logged_9var_processed,by='feature')
var9_all = full_join(var9_all,nonlogged_9var_processed,by='feature')
var9_all[is.na(var9_all)] = 0

var_all = full_join(var3_all,var6_all,by='feature')
var_all = full_join(var_all,var9_all,by='feature') %>% select(-feature)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

outmap = var_all %>% cor %>% pheatmap
save_pheatmap_pdf(outmap, paste(args[[1]],'_janus_heatmap.pdf',sep=''))

