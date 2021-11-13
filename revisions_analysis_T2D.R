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
cohort2_abundances = voe_out$original_data$dependent_variables[[2]] 
cohort3_abundances = voe_out$original_data$dependent_variables[[3]] 

cohort1_metadata = voe_out$original_data$independent_variables[[1]]
cohort2_metadata = voe_out$original_data$independent_variables[[2]]
cohort3_metadata = voe_out$original_data$independent_variables[[3]]

# get subsampled columns to consider per cohort
columns_to_consider_cohort1 = cohort1_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort2 = cohort2_metadata %>% select(-sampleID,-study_condition) %>% colnames
columns_to_consider_cohort3 = cohort3_metadata %>% select(-sampleID,-study_condition) %>% colnames

col3co1 = sample(columns_to_consider_cohort1,3)
col3co2 = sample(columns_to_consider_cohort2,3)
col3co3 = sample(columns_to_consider_cohort3,3)

col6co1 = sample(columns_to_consider_cohort1,6)
col6co2 = sample(columns_to_consider_cohort2,6)
col6co3 = sample(columns_to_consider_cohort3,6)

col9co1 = sample(columns_to_consider_cohort1,9)
col9co2 = sample(columns_to_consider_cohort2,9)
col9co3 = sample(columns_to_consider_cohort3,9)

# 3 varables in total
col3mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col3co1))
col3mdat2 = cohort2_metadata %>% select(sampleID,study_condition,all_of(col3co2))
col3mdat3 = cohort3_metadata %>% select(sampleID,study_condition,all_of(col3co3))

# 6 variables in total
col6mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col6co1))
col6mdat2 = cohort2_metadata %>% select(sampleID,study_condition,all_of(col6co2))
col6mdat3 = cohort3_metadata %>% select(sampleID,study_condition,all_of(col6co3))

# 9 variables in total
col9mdat1 = cohort1_metadata %>% select(sampleID,study_condition,all_of(col9co1))
col9mdat2 = cohort2_metadata %>% select(sampleID,study_condition,all_of(col9co2))
col9mdat3 = cohort3_metadata %>% select(sampleID,study_condition,all_of(col9co3))

# meta analysis datasets for each var number
col3_meta_analysis = list(col3mdat1,col3mdat2,col3mdat3)
col6_meta_analysis = list(col6mdat1,col6mdat2,col6mdat3)
col9_meta_analysis = list(col9mdat1,col9mdat2,col9mdat3)

### logged voe
voe_df_3var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat2_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat2, dependent_variables = cohort2_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat3_logged = quantvoe::full_voe_pipeline(independent_variables = col3mdat3, dependent_variables = cohort3_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_meta = quantvoe::full_voe_pipeline(independent_variables = col3_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_3var_mdat1_logged,voe_df_3var_mdat2_logged,voe_df_3var_mdat3_logged,voe_df_3var_meta_logged),'logged_3var_voe.rds')

voe_df_6var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat2_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat2, dependent_variables = cohort2_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat3_logged = quantvoe::full_voe_pipeline(independent_variables = col6mdat3, dependent_variables = cohort3_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_meta = quantvoe::full_voe_pipeline(independent_variables = col6_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_6var_mdat1_logged,voe_df_6var_mdat2_logged,voe_df_6var_mdat3_logged,voe_df_6var_meta_logged),'logged_6var_voe.rds')

voe_df_9var_mdat1_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat2_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat2, dependent_variables = cohort2_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat3_logged = quantvoe::full_voe_pipeline(independent_variables = col9mdat3, dependent_variables = cohort3_abundances, primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_meta = quantvoe::full_voe_pipeline(independent_variables = col9_meta_analysis, dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_9var_mdat1_logged,voe_df_9var_mdat2_logged,voe_df_9var_mdat3_logged,voe_df_9var_meta_logged),'logged_9var_voe.rds')

### non-logged voe
cohort1_abundances_nonlogged = voe_out$original_data$dependent_variables[[1]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort2_abundances_nonlogged = voe_out$original_data$dependent_variables[[2]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))
cohort3_abundances_nonlogged = voe_out$original_data$dependent_variables[[3]] %>% mutate_if(is.numeric,function(x) round(exp(x) - 0.0000000100,20))

voe_df_3var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat2_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat2, dependent_variables = cohort2_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat3_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3mdat3, dependent_variables = cohort3_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_meta_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col3_meta_analysis, dependent_variables = list(cohort1_abundances_nonlogged,cohort2_abundances_nonlogged,cohort3_abundances_nonlogged), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_3var_mdat1_nonlogged,voe_df_3var_mdat2_nonlogged,voe_df_3var_mdat3_nonlogged,voe_df_3var_meta_nonlogged),'nonlogged_3var_voe.rds')

voe_df_6var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat2_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat2, dependent_variables = cohort2_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat3_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6mdat3, dependent_variables = cohort3_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_meta_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col6_meta_analysis, dependent_variables = list(cohort1_abundances_nonlogged,cohort2_abundances_nonlogged,cohort3_abundances_nonlogged), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_6var_mdat1_nonlogged,voe_df_6var_mdat2_nonlogged,voe_df_6var_mdat3_nonlogged,voe_df_6var_meta_nonlogged),'nonlogged_6var_voe.rds')

voe_df_9var_mdat1_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat2_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat2, dependent_variables = cohort2_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat3_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9mdat3, dependent_variables = cohort3_abundances_nonlogged , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_meta_nonlogged = quantvoe::full_voe_pipeline(independent_variables = col9_meta_analysis, dependent_variables = list(cohort1_abundances_nonlogged,cohort2_abundances_nonlogged,cohort3_abundances_nonlogged), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_9var_mdat1_nonlogged,voe_df_9var_mdat2_nonlogged,voe_df_9var_mdat3_nonlogged,voe_df_9var_meta_nonlogged),'nonlogged_9var_voe.rds')

#### CLR transformation voe
cohort1_abundances_nonlogged[,2:ncol(cohort1_abundances_nonlogged)] = clr(cohort1_abundances_nonlogged[,2:ncol(cohort1_abundances_nonlogged)])
cohort1_abundances_clr = cohort1_abundances_nonlogged
cohort2_abundances_nonlogged[,2:ncol(cohort2_abundances_nonlogged)] = clr(cohort2_abundances_nonlogged[,2:ncol(cohort2_abundances_nonlogged)])
cohort2_abundances_clr = cohort2_abundances_nonlogged
cohort3_abundances_nonlogged[,2:ncol(cohort3_abundances_nonlogged)] = clr(cohort3_abundances_nonlogged[,2:ncol(cohort3_abundances_nonlogged)])
cohort3_abundances_clr = cohort3_abundances_nonlogged

voe_df_3var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat1, dependent_variables = cohort1_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat2_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat2, dependent_variables = cohort2_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_mdat3_clr = quantvoe::full_voe_pipeline(independent_variables = col3mdat3, dependent_variables = cohort3_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_3var_meta_clr = quantvoe::full_voe_pipeline(independent_variables = col3_meta_analysis, dependent_variables = list(cohort1_abundances_clr,cohort2_abundances_clr,cohort3_abundances_clr), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_3var_mdat1_clr,voe_df_3var_mdat2_clr,voe_df_3var_mdat3_clr,voe_df_3var_meta_clr),'clr_3var_voe.rds')

voe_df_6var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat1, dependent_variables = cohort1_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat2_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat2, dependent_variables = cohort2_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_mdat3_clr = quantvoe::full_voe_pipeline(independent_variables = col6mdat3, dependent_variables = cohort3_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_6var_meta_clr = quantvoe::full_voe_pipeline(independent_variables = col6_meta_analysis, dependent_variables = list(cohort1_abundances_clr,cohort2_abundances_clr,cohort3_abundances_clr), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_6var_mdat1_clr,voe_df_6var_mdat2_clr,voe_df_6var_mdat3_clr,voe_df_6var_meta_clr),'clr_6var_voe.rds')

voe_df_9var_mdat1_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat1, dependent_variables = cohort1_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat2_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat2, dependent_variables = cohort2_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_mdat3_clr = quantvoe::full_voe_pipeline(independent_variables = col9mdat3, dependent_variables = cohort3_abundances_clr , primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1)
voe_df_9var_meta_clr = quantvoe::full_voe_pipeline(independent_variables = col9_meta_analysis, dependent_variables = list(cohort1_abundances_clr,cohort2_abundances_clr,cohort3_abundances_clr), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

saveRDS(list(voe_df_9var_mdat1_clr,voe_df_9var_mdat2_clr,voe_df_9var_mdat3_clr,voe_df_9var_meta_clr),'clr_9var_voe.rds')

#### subsampling number of vibrations analysis
voe_df_all_var_meta10 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=10, fdr_cutoff = 1, meta_analysis=TRUE)
voe_df_all_var_meta100 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=100, fdr_cutoff = 1, meta_analysis=TRUE)
voe_df_all_var_meta500 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=500, fdr_cutoff = 1, meta_analysis=TRUE)
voe_df_all_var_meta1000 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=1000, fdr_cutoff = 1, meta_analysis=TRUE)
voe_df_all_var_meta2500 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=2500, fdr_cutoff = 1, meta_analysis=TRUE)
voe_df_all_var_meta5000 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=5000, fdr_cutoff = 1, meta_analysis=TRUE)
voe_df_all_var_meta7500 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=7500, fdr_cutoff = 1, meta_analysis=TRUE)
voe_df_all_var_meta10000 = quantvoe::full_voe_pipeline(independent_variables = list(cohort1_metadata,cohort2_metadata,cohort3_metadata), dependent_variables = list(cohort1_abundances,cohort2_abundances,cohort3_abundances), primary_variable = 'study_condition', max_vibration_num=10000, fdr_cutoff = 1, meta_analysis=TRUE)

data = list(voe_df_all_var_meta10,voe_df_all_var_meta100,voe_df_all_var_meta500,voe_df_all_var_meta1000,voe_df_all_var_meta2500,voe_df_all_var_meta5000,voe_df_all_var_meta7500,voe_df_all_var_meta10000)
data = map(data, function(x) x['vibration_output']$vibration_output$confounder_analysis$summarized_output)

saveRDS(data,'vibration_subsampling_analysis.rds')











