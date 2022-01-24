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

voe_out = readRDS('T2D_90.rds')

# load abundances and reset to non logged
cohort1_abundances = voe_out$original_data$dependent_variables[[1]] 

cohort1_metadata = voe_out$original_data$independent_variables[[1]]

# get subsampled columns to consider per cohort
columns_to_consider_cohort1 = cohort1_metadata %>% select(-sampleID,-study_condition) %>% colnames

### ADDING ADDITION SAMPLING FOR COHORT WITH THE MOST VARIABLES PER R2's SECOND SET OF REQUESTS
col12co1 = sample(columns_to_consider_cohort1,12)
col15co1 = sample(columns_to_consider_cohort1,15)
col18co1 = sample(columns_to_consider_cohort1,18)
col21co1 = sample(columns_to_consider_cohort1,21)
col24co1 = sample(columns_to_consider_cohort1,24)

### ADDING ADDITION SAMPLING FOR COHORT WITH THE MOST VARIABLES PER R2's SECOND SET OF REQUESTS
col12mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col12co1))
col15mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col15co1))
col18mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col18co1))
col21mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col21co1))
col24mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col24co1))

### ADDING ADDITION SAMPLING FOR COHORT WITH THE MOST VARIABLES PER R2's SECOND SET OF REQUESTS
voe_df_12var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col12mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_15var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col15mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_18var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col18mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_21var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col21mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_24var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col24p;;;;mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_12var_mdat1_logged,voe_df_15var_mdat1_logged,voe_df_18var_mdat1_logged,voe_df_21var_mdat1_logged,voe_df_24var_mdat1_logged),'logged_cohort1_12_to_24_var_voe_T2D.rds')

### ADDING ADDITION SAMPLING FOR COHORT WITH THE MOST VARIABLES PER R2's SECOND SET OF REQUESTS
voe_df_12var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col12mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_15var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col15mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_18var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col18mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_21var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col21mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_24var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col24mdat1, dependent_variables = cohort1_abundances_nonlogged, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_12var_mdat1_nonlogged,voe_df_15var_mdat1_nonlogged,voe_df_18var_mdat1_nonlogged,voe_df_21var_mdat1_nonlogged,voe_df_24var_mdat1_nonlogged),'nonlogged_cohort1_12_to_24_var_voe_T2D.rds')

### ADDING ADDITION SAMPLING FOR COHORT WITH THE MOST VARIABLES PER R2's SECOND SET OF REQUESTS
voe_df_12var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col12mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_15var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col15mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_18var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col18mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_21var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col21mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_24var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col24mdat1, dependent_variables = cohort1_abundances_clr, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)

saveRDS(list(voe_df_12var_mdat1_clr,voe_df_15var_mdat1_clr,voe_df_18var_mdat1_clr,voe_df_21var_mdat1_clr,voe_df_24var_mdat1_clr),'clr_cohort1_12_to_24_var_voe_T2D.rds')









