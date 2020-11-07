#preprocess curated metagenomic data for pipeline

library(tidyverse)

###average all abundance data and create overall file

metadata = readRDS('cmd_all_metadata_2.rds')
abundance_data = t(readRDS('metaphlan_abundance_cmd.rds')) %>% data.frame %>% rownames_to_column('sampleID') %>% mutate(sampleID = strsplit(sampleID,':') %>% map_chr(2))

sample_subject_dataset_mapping  = metadata %>% select(sampleID,subjectID,dataset_name)
abundance_data = inner_join(sample_subject_dataset_mapping,abundance_data) %>% select(-sampleID)

average_feature = function(feature_name,data){
	data = data %>% select(dataset_name,subjectID,feature_name) %>% group_by(dataset_name,subjectID) %>% summarise(.groups='keep',taxa = mean(get(feature_name))) %>% distinct
	colnames(data)[3]=feature_name
	data = data %>% ungroup %>% arrange(subjectID,dataset_name) %>% select(-dataset_name)
	return(data)
}

output = list()
microbes=colnames(abundance_data %>% select(-dataset_name,-subjectID))
for(taxa in microbes){
	print(taxa)
	out = average_feature(taxa,abundance_data) 
	subs = out %>% select(subjectID)
	output[[taxa]] = out %>% select(-subjectID)
}

averaged_data = bind_cols(subs,output)

saveRDS(averaged_data,'curated_metagenomic_data_averaged_redo.rds')

###get metadata dataframes by dataset
datasets = unique(sample_subject_dataset_mapping$dataset_name)


to_remove_initial <- c("PMID", "NCBI_accession", "number_bases", "curator", "sequencing_platform", "DNA_extraction_kit", "body_site", "disease", "non_westernized", "location","minimum_read_length", "median_read_length", "antibiotics_current_use") 

for(p in c('CRC','cirrhosis','IBD','T2D','ACVD','T1D')){
	system(paste('rm -rf',p))
	system(paste('mkdir',p))
}

#for each dataset, subset
for(d in datasets){
	metadata_sub = metadata %>% filter(dataset_name == d)  %>% select(-sampleID)
	#for each dataset, get the phenotypes (cases and controls)
	phenotypes = metadata_sub %>% filter(study_condition != 'control') %>% select(study_condition) %>% unique %>% unlist %>% unname
	phenotypes = intersect(c('CRC','cirrhosis','IBD','T2D','ACVD','T1D'),phenotypes)
	#if in the 6 phenotypes of interest
	if(length(phenotypes)!=0){
		for(p in phenotypes){
			print(p)
			metadata_sub = metadata_sub %>% filter(study_condition == 'control' | study_condition == p) %>% select(-all_of(to_remove_initial))
			#get number of subjects and number of samples
			subnum = length(unique(metadata_sub$subjectID))
			samnum = length(unique(metadata_sub$sampleID))
			#if subjects<samples need to collapse metadata along subjects
			if(subnum!=samnum){
				subs = metadata_sub %>% select(subjectID) 
				#get numeric columns and avereage
				numeric = metadata_sub %>% select_if(is.numeric) 
				numeric = bind_cols(subs,numeric) %>% group_by(subjectID) %>% summarise_each(mean)
				#get character columns and collapse
				character = metadata_sub %>% select_if(is.character)%>% distinct
	            toremove=list()
	            for(s in subs %>% unique %>% unlist %>% unname){
	              character_sub = character %>% filter(subjectID==s)
	              temp=apply(character_sub, 2, function(x) length(unique(x)))
	              temp=names(temp[temp>1])
	              if(length(temp)>0){
	                toremove[[s]]=temp
	         	   }
	        	}
	            toremove=unique(unname(unlist(toremove)))
	            if(length(toremove)>0){
	              character=character %>% select(-all_of(toremove)) %>% distinct()
	            }
				metadata_sub = inner_join(numeric,character)
			}
			metadata_sub =  metadata_sub[,colSums(is.na(metadata_sub))/nrow(metadata_sub)<.1]  %>% select_if(~n_distinct(.) > 1)
			#retrieve relevant rows from abundance data and log them
			abundance_sub = averaged_data %>% filter(subjectID %in% metadata_sub$subjectID) %>% distinct()
			ab_subs = abundance_sub %>% select(subjectID)
			abs = abundance_sub %>% select(-subjectID)
			abs = log(0.00000001+abs)
			abundance_sub = bind_cols(ab_subs,abs)
			saveRDS(abundance_sub,paste(p,'/',d,'_',p,'_','abundance_repr.rds',sep=''))
			saveRDS(metadata_sub,paste(p,'/',d,'_',p,'_','metadata_repr.rds',sep=''))
		}
	}
}






