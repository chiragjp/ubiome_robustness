# compare impact on number of vibrations on confounder analysis

library(tidyverse)
library(pheatmap)

setwd('~/Dropbox (HMS)/RagGroup Team Folder/Braden Tierney/vibration_of_effects_microbiome/repr/revisions/')

files = list.files()

files = files[grep('subsampling',files)]

output = list()
for(f in files){
  data = readRDS(f)
  name = gsub('vibration_subsampling_analysis_','',f)
  name = gsub('.rds','',name)
  data = data[[1]] %>% mutate(number_of_vibrations = as.numeric(name))
  output[[f]] = data
}

output_merged = bind_rows(output) %>% select(term,estimate,number_of_vibrations) %>% filter(term != '(Intercept)',term != 'sd__(Intercept)', term!= 'sd__Observation') 
forheat1 = dcast(output_merged, term ~ number_of_vibrations,value.var = "estimate") %>% column_to_rownames("term") 
forheat2 = forheat1 %>% cor

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

outmap1 = pheatmap(forheat1,cluster_cols = F,cluster_rows = F)

save_pheatmap_pdf(outmap1, 't2d_association_size_heatmap.pdf')

outmap2 = pheatmap(forheat2,cluster_cols = F,cluster_rows = F)

save_pheatmap_pdf(outmap2, 't2d_vibration_heatmap.pdf')



