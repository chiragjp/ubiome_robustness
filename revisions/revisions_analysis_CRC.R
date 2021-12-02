# revision analysis, subsampling T2D

### analyses
# Try different logging strategies -- code written
# Different numbers of variables -- code written 
# Different numbers of vibrations -- code written

### figures to generate (in assocated RMD)

# upset R plots, number of features significant etc for each:
	# logging method
	# number of variables
	# cohort/MA strat

# grids of voe plots for all features across different var nums + different cohorts/ma
	# logging method on columns, cohorts on rows, different figure for variable number

# correlation between number of vibrations and direction of association

library(tidyverse)
library(compositions)
library(quantvoe)

voe_out = readRDS('CRC_90.rds')

# load abundances and reset to non logged
cohort1_abundances = voe_out$original_data$dependent_variables[[1]] 
cohort2_abundances = voe_out$original_data$dependent_variables[[2]] 
cohort3_abundances = voe_out$original_data$dependent_variables[[3]] 
cohort4_abundances = voe_out$original_data$dependent_variables[[4]] 
cohort5_abundances = voe_out$original_data$dependent_variables[[5]] 
cohort6_abundances = voe_out$original_data$dependent_variables[[6]] 
cohort7_abundances = voe_out$original_data$dependent_variables[[7]] 

cohort1_metadata = voe_out$original_data$independent_variables[[1]]
cohort2_metadata = voe_out$original_data$independent_variables[[2]]
cohort3_metadata = voe_out$original_data$independent_variables[[3]]
cohort4_metadata = voe_out$original_data$independent_variables[[4]]
cohort5_metadata = voe_out$original_data$independent_variables[[5]]
cohort6_metadata = voe_out$original_data$independent_variables[[6]]
cohort7_metadata = voe_out$original_data$independent_variables[[7]]

# get subsampled columns to consider per cohort
columns_to_consider_cohort1 = cohort1_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort2 = cohort2_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort3 = cohort3_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort4 = cohort4_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort5 = cohort5_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort6 = cohort6_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort7 = cohort7_metadata %>% select(-sampleID,-study_condition) %>% colnames

col3co1 = sample(columns_to_consider_cohort1,3)
col3co2 = sample(columns_to_consider_cohort2,3)
col3co3 = sample(columns_to_consider_cohort3,3)
col3co4 = sample(columns_to_consider_cohort4,3)
col3co5 = sample(columns_to_consider_cohort5,3)
col3co6 = sample(columns_to_consider_cohort6,3)
col3co7 = sample(columns_to_consider_cohort7,3)

col6co1 = sample(columns_to_consider_cohort1,6,replace=TRUE) %>% unique
col6co2 = sample(columns_to_consider_cohort2,6,replace=TRUE) %>% unique
col6co3 = sample(columns_to_consider_cohort3,6,replace=TRUE) %>% unique
col6co4 = sample(columns_to_consider_cohort4,6,replace=TRUE) %>% unique
col6co5 = sample(columns_to_consider_cohort5,6,replace=TRUE) %>% unique
col6co6 = sample(columns_to_consider_cohort6,6,replace=TRUE) %>% unique
col6co7 = sample(columns_to_consider_cohort7,6,replace=TRUE) %>% unique

col9co1 = sample(columns_to_consider_cohort1,9,replace=TRUE) %>% unique
col9co2 = sample(columns_to_consider_cohort2,9,replace=TRUE) %>% unique
col9co3 = sample(columns_to_consider_cohort3,9,replace=TRUE) %>% unique
col9co4 = sample(columns_to_consider_cohort4,9,replace=TRUE) %>% unique
col9co5 = sample(columns_to_consider_cohort5,9,replace=TRUE) %>% unique
col9co6 = sample(columns_to_consider_cohort6,9,replace=TRUE) %>% unique
col9co7 = sample(columns_to_consider_cohort7,9,replace=TRUE) %>% unique

# 3 varables in total
col3mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col3co1))
col3mdat2 = cohort2_metadata %>% select(sampleID,study_condition,all_of(col3co2))
col3mdat3 = cohort3_metadata %>% select(sampleID,study_condition,all_of(col3co3))
col3mdat4 = cohort4_metadata %>% select(sampleID,study_condition,all_of(col3co4))
col3mdat5 = cohort5_metadata %>% select(sampleID,study_condition,all_of(col3co5))
col3mdat6 = cohort6_metadata %>% select(sampleID,study_condition,all_of(col3co6))
col3mdat7 = cohort7_metadata %>% select(sampleID,study_condition,all_of(col3co7))

# 6 variables in total
col6mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col6co1))
col6mdat2 = cohort2_metadata %>% select(sampleID,study_condition,all_of(col6co2))
col6mdat3 = cohort3_metadata %>% select(sampleID,study_condition,all_of(col6co3))
col6mdat4 = cohort4_metadata %>% select(sampleID,study_condition,all_of(col6co4))
col6mdat5 = cohort5_metadata %>% select(sampleID,study_condition,all_of(col6co5))
col6mdat6 = cohort6_metadata %>% select(sampleID,study_condition,all_of(col6co6))
col6mdat7 = cohort7_metadata %>% select(sampleID,study_condition,all_of(col6co7))

# 9 variables in total
col9mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col9co1))
col9mdat2 = cohort2_metadata %>% select(sampleID,study_condition,all_of(col9co2))
col9mdat3 = cohort3_metadata %>% select(sampleID,study_condition,all_of(col9co3))
col9mdat4 = cohort4_metadata %>% select(sampleID,study_condition,all_of(col9co4))
col9mdat5 = cohort5_metadata %>% select(sampleID,study_condition,all_of(col9co5))
col9mdat6 = cohort6_metadata %>% select(sampleID,study_condition,all_of(col9co6))
col9mdat7 = cohort7_metadata %>% select(sampleID,study_condition,all_of(col9co7))

# meta analysis datasets for each var number
col3_meta_analysis = list(col3mdat1,col3mdat2,col3mdat3,col3mdat4,col3mdat5,col3mdat6,col3mdat7)
col6_meta_analysis = list(col6mdat1,col6mdat2,col6mdat3,col6mdat4,col6mdat5,col6mdat6,col6mdat7)
col9_meta_analysis = list(col9mdat1,col9mdat2,col9mdat3,col9mdat4,col9mdat5,col9mdat6,col9mdat7)

### logged voe
voe_df_3var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat2_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat2, dependent_variables = cohort2_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat3_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat3, dependent_variables = cohort3_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat4_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat4, dependent_variables = cohort4_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat5_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat5, dependent_variables = cohort5_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat6_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat6, dependent_variables = cohort6_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat7_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat7, dependent_variables = cohort7_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_meta_logged = quantvoe::full_voe_pipeline(independent_variables = col3_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_3var_mdat1_logged,voe_df_3var_mdat2_logged,voe_df_3var_mdat3_logged,voe_df_3var_mdat4_logged,voe_df_3var_mdat5_logged,voe_df_3var_mdat6_logged,voe_df_3var_mdat7_logged,voe_df_3var_meta_logged),'logged_3var_voe_CRC.rds')

voe_df_6var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat2_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat2, dependent_variables = cohort2_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat3_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat3, dependent_variables = cohort3_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat4_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat4, dependent_variables = cohort4_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat5_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat5, dependent_variables = cohort5_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat6_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat6, dependent_variables = cohort6_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat7_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat7, dependent_variables = cohort7_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_meta_logged = quantvoe::full_voe_pipeline(independent_variables = col6_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_6var_mdat1_logged,voe_df_6var_mdat2_logged,voe_df_6var_mdat3_logged,voe_df_6var_mdat4_logged,voe_df_6var_mdat5_logged,voe_df_6var_mdat6_logged,voe_df_6var_mdat7_logged,voe_df_6var_meta_logged),'logged_6var_voe_CRC.rds')

voe_df_9var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat2_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat2, dependent_variables = cohort2_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat3_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat3, dependent_variables = cohort3_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat4_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat4, dependent_variables = cohort4_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat5_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat5, dependent_variables = cohort5_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat6_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat6, dependent_variables = cohort6_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat7_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat7, dependent_variables = cohort7_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_meta_logged = quantvoe::full_voe_pipeline(independent_variables = col9_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_9var_mdat1_logged,voe_df_9var_mdat2_logged,voe_df_9var_mdat3_logged,voe_df_9var_mdat4_logged,voe_df_9var_mdat5_logged,voe_df_9var_mdat6_logged,voe_df_9var_mdat7_logged,voe_df_9var_meta_logged),'logged_9var_voe_CRC.rds')

### non-logged voe
cohort1_abundances_nonlogged = voe_out$original_data$dependent_variables[[1]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort2_abundances_nonlogged = voe_out$original_data$dependent_variables[[2]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort3_abundances_nonlogged = voe_out$original_data$dependent_variables[[3]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort4_abundances_nonlogged = voe_out$original_data$dependent_variables[[4]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort5_abundances_nonlogged = voe_out$original_data$dependent_variables[[5]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort6_abundances_nonlogged = voe_out$original_data$dependent_variables[[6]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort7_abundances_nonlogged = voe_out$original_data$dependent_variables[[7]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))

voe_df_3var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat2_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat2, dependent_variables = cohort2_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat3_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat3, dependent_variables = cohort3_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat4_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat4, dependent_variables = cohort4_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat5_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat5, dependent_variables = cohort5_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat6_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat6, dependent_variables = cohort6_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat7_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat7, dependent_variables = cohort7_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_meta_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_3var_mdat1_nonlogged,voe_df_3var_mdat2_nonlogged,voe_df_3var_mdat3_nonlogged,voe_df_3var_mdat4_nonlogged,voe_df_3var_mdat5_nonlogged,voe_df_3var_mdat6_nonlogged,voe_df_3var_mdat7_nonlogged,voe_df_3var_meta_nonlogged),'nonlogged_3var_voe_CRC.rds')

voe_df_6var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat2_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat2, dependent_variables = cohort2_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat3_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat3, dependent_variables = cohort3_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat4_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat4, dependent_variables = cohort4_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat5_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat5, dependent_variables = cohort5_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat6_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat6, dependent_variables = cohort6_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat7_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat7, dependent_variables = cohort7_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_meta_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_6var_mdat1_nonlogged,voe_df_6var_mdat2_nonlogged,voe_df_6var_mdat3_nonlogged,voe_df_6var_mdat4_nonlogged,voe_df_6var_mdat5_nonlogged,voe_df_6var_mdat6_nonlogged,voe_df_6var_mdat7_nonlogged,voe_df_6var_meta_nonlogged),'nonlogged_6var_voe_CRC.rds')

voe_df_9var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat2_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat2, dependent_variables = cohort2_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat3_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat3, dependent_variables = cohort3_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat4_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat4, dependent_variables = cohort4_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat5_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat5, dependent_variables = cohort5_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat6_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat6, dependent_variables = cohort6_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat7_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat7, dependent_variables = cohort7_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_meta_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_9var_mdat1_nonlogged,voe_df_9var_mdat2_nonlogged,voe_df_9var_mdat3_nonlogged,voe_df_9var_mdat4_nonlogged,voe_df_9var_mdat5_nonlogged,voe_df_9var_mdat6_nonlogged,voe_df_9var_mdat7_nonlogged,voe_df_9var_meta_nonlogged),'nonlogged_9var_voe_CRC.rds')

#### CLR transformation voe
cohort1_abundances_nonlogged[,2:ncol(cohort1_abundances_nonlogged)] = clr(cohort1_abundances_nonlogged[,2:ncol(cohort1_abundances_nonlogged)])
cohort1_abundances_clr = cohort1_abundances_nonlogged
cohort2_abundances_nonlogged[,2:ncol(cohort2_abundances_nonlogged)] = clr(cohort2_abundances_nonlogged[,2:ncol(cohort2_abundances_nonlogged)])
cohort2_abundances_clr = cohort2_abundances_nonlogged
cohort3_abundances_nonlogged[,2:ncol(cohort3_abundances_nonlogged)] = clr(cohort3_abundances_nonlogged[,2:ncol(cohort3_abundances_nonlogged)])
cohort3_abundances_clr = cohort3_abundances_nonlogged
cohort4_abundances_nonlogged[,2:ncol(cohort4_abundances_nonlogged)] = clr(cohort4_abundances_nonlogged[,2:ncol(cohort4_abundances_nonlogged)])
cohort4_abundances_clr = cohort4_abundances_nonlogged
cohort5_abundances_nonlogged[,2:ncol(cohort5_abundances_nonlogged)] = clr(cohort5_abundances_nonlogged[,2:ncol(cohort5_abundances_nonlogged)])
cohort5_abundances_clr = cohort5_abundances_nonlogged
cohort6_abundances_nonlogged[,2:ncol(cohort6_abundances_nonlogged)] = clr(cohort6_abundances_nonlogged[,2:ncol(cohort6_abundances_nonlogged)])
cohort6_abundances_clr = cohort6_abundances_nonlogged
cohort7_abundances_nonlogged[,2:ncol(cohort7_abundances_nonlogged)] = clr(cohort7_abundances_nonlogged[,2:ncol(cohort7_abundances_nonlogged)])
cohort7_abundances_clr = cohort7_abundances_nonlogged

voe_df_3var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat2_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat2, dependent_variables = cohort2_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat3_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat3, dependent_variables = cohort3_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat4_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat4, dependent_variables = cohort4_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat5_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat5, dependent_variables = cohort5_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat6_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat6, dependent_variables = cohort6_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat7_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat7, dependent_variables = cohort7_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_meta_clr = quantvoe::full_voe_pipeline(independent_variables = col3_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_3var_mdat1_clr,voe_df_3var_mdat2_clr,voe_df_3var_mdat3_clr,voe_df_3var_mdat4_clr,voe_df_3var_mdat5_clr,voe_df_3var_mdat6_clr,voe_df_3var_mdat7_clr,voe_df_3var_meta_clr),'clr_3var_voe_CRC.rds')

voe_df_6var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat2_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat2, dependent_variables = cohort2_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat3_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat3, dependent_variables = cohort3_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat4_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat4, dependent_variables = cohort4_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat5_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat5, dependent_variables = cohort5_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat6_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat6, dependent_variables = cohort6_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat7_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat7, dependent_variables = cohort7_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_meta_clr = quantvoe::full_voe_pipeline(independent_variables = col6_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_6var_mdat1_clr,voe_df_6var_mdat2_clr,voe_df_6var_mdat3_clr,voe_df_6var_mdat4_clr,voe_df_6var_mdat5_clr,voe_df_6var_mdat6_clr,voe_df_6var_mdat7_clr,voe_df_6var_meta_clr),'clr_6var_voe_CRC.rds')

voe_df_9var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat2_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat2, dependent_variables = cohort2_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat3_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat3, dependent_variables = cohort3_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat4_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat4, dependent_variables = cohort4_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat5_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat5, dependent_variables = cohort5_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat6_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat6, dependent_variables = cohort6_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat7_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat7, dependent_variables = cohort7_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_meta_clr = quantvoe::full_voe_pipeline(independent_variables = col9_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances,cohort4_abundances,cohort5_abundances,cohort6_abundances,cohort7_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)
saveRDS(list(voe_df_9var_mdat1_clr,voe_df_9var_mdat2_clr,voe_df_9var_mdat3_clr,voe_df_9var_mdat4_clr,voe_df_9var_mdat5_clr,voe_df_9var_mdat6_clr,voe_df_9var_mdat7_clr,voe_df_9var_meta_clr),'clr_9var_voe_CRC.rds')












