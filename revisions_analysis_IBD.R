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

voe_out = readRDS('IBD_90.rds')

# load abundances and reset to non logged
cohort1_abundances = voe_out$original_data$dependent_variables[[1]] 

cohort1_metadata = voe_out$original_data$independent_variables[[1]]

# get subsampled columns to consider per cohort
columns_to_consider_cohort1 = cohort1_metadata %>% select(-sampleID,-study_condition) %>% colnames

col3co1 = sample(columns_to_consider_cohort1,2)

col6co1 = sample(columns_to_consider_cohort1,3)

col9co1 = sample(columns_to_consider_cohort1,4)

# 3 varables in total
col3mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col3co1))

# 6 variables in total
col6mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col6co1))

# 9 variables in total
col9mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col9co1))


### logged voe
voe_df_3var_mdat1 = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_3var_mdat1_logged),'logged_2var_voe_IBD.rds')

voe_df_6var_mdat1 = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_6var_mdat1_logged),'logged_3var_voe_IBD.rds')

voe_df_9var_mdat1 = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_9var_mdat1_logged),'logged_4var_voe_IBD.rds')

### non-logged voe
cohort1_abundances_nonlogged = voe_out$original_data$dependent_variables[[1]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))

voe_df_3var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_3var_mdat1_nonlogged),'nonlogged_2var_voe_IBD.rds')

voe_df_6var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_6var_mdat1_nonlogged),'nonlogged_3var_voe_IBD.rds')

voe_df_9var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_9var_mdat1_nonlogged),'nonlogged_4var_voe_IBD.rds')

#### CLR transformation voe
cohort1_abundances_nonlogged[,2:ncol(cohort1_abundances_nonlogged)] = clr(cohort1_abundances_nonlogged[,2:ncol(cohort1_abundances_nonlogged)])
cohort1_abundances_clr = cohort1_abundances_nonlogged

voe_df_3var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_3var_mdat1_clr),'clr_2var_voe_IBD.rds')

voe_df_6var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_6var_mdat1_clr),'clr_3var_voe_IBD.rds')

voe_df_9var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_9var_mdat1_clr),'clr_4var_voe_IBD.rds')








