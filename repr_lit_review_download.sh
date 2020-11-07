#!/bin/bash

esearch -db pubmed --query "\"type 2 diabetes\" AND \"gut microbiome\"" | efilter -pub journal | efetch -format csv > all_t2d_pubs2.csv
esearch -db pubmed --query "\"type 1 diabetes\" AND \"gut microbiome\"" | efilter -pub journal | efetch -format csv> all_t1d_pubs.csv
esearch -db pubmed --query "\"colorectal cancer\" OR CRC AND \"gut microbiome\"" | efilter -pub journal | efetch -format csv> all_crc_pubs.csv
esearch -db pubmed --query "ACVD OR atherosclerosis OR \"atherosclerotic cardiovascular disease\" OR \"cardiovascular disease\" AND \"gut microbiome\" " | efilter -pub journal | efetch -format csv> all_acvd_pubs.csv
esearch -db pubmed --query "cirrhosis AND \"gut microbiome\"" | efilter -pub journal | efetch -format csv> all_cirrhosis_pubs.csv
esearch -db pubmed --query "IBD OR \"inflammatory bowel disease\" OR crohn's disease OR \"ulcerative colitis\" AND \"gut microbiome\"" | efilter -pub journal | efetch -format csv> all_ibd_pubs.csv
