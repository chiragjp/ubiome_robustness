library(tidyverse)
library(ComplexUpset)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

files = read.csv(args[[1]],header=F) %>% unlist %>% unname

clr_3var = readRDS(files[[1]])
clr_3var[[4]] = clr_3var[[4]][names(clr_3var[[4]])!='initial_association_output']
clr_6var = readRDS(files[[2]])
clr_6var[[4]] = clr_6var[[4]][names(clr_6var[[4]])!='initial_association_output']
clr_9var = readRDS(files[[3]])
clr_9var[[4]] = clr_9var[[4]][names(clr_9var[[4]])!='initial_association_output']

logged_3var = readRDS(files[[4]])
logged_3var[[4]] = logged_3var[[4]][names(logged_3var[[4]])!='initial_association_output']
logged_6var = readRDS(files[[5]])
logged_6var[[4]] = logged_6var[[4]][names(logged_6var[[4]])!='initial_association_output']
logged_9var = readRDS(files[[6]])
logged_9var[[4]] = logged_9var[[4]][names(logged_9var[[4]])!='initial_association_output']

nonlogged_3var = readRDS(files[[7]])
nonlogged_3var[[4]] = nonlogged_3var[[4]][names(nonlogged_3var[[4]])!='initial_association_output']
nonlogged_6var = readRDS(files[[8]])
nonlogged_6var[[4]] = nonlogged_6var[[4]][names(nonlogged_6var[[4]])!='initial_association_output']
nonlogged_9var = readRDS(files[[9]])
nonlogged_9var[[4]] = nonlogged_9var[[4]][names(nonlogged_9var[[4]])!='initial_association_output']
