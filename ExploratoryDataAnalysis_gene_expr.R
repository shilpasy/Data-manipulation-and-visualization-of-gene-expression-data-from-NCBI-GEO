#Data wrangling - gene expression data
library(dplyr)
library(tidyverse)
library(GEOquery)

#setwd('D:/shilpa/GSE')
#Data was obtained from:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947

data <- read.csv(file = "GSE183947_fpkm.csv")

gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

#Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# select, mutate (edit), rename ------------
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>% #global substitution equivalent to replace in python
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

#head(data)

# reshaping data - basically converting the matrix format into columns (long) so that you can join it with above metadata)
#combining with metadata will help in our analysis

data.long <- data %>%
  rename(gene = X) %>% 
  gather(key = 'samples', value = 'FPKM', -gene)


data.long <- data.long %>%
  left_join(., metadata.modified, by = c("samples" = "description")) 

write.csv(data.long, 'gene_samples_FPKM.csv', row.names = FALSE)  #this file can be used to exploratory data analysis

# explore data ## Summary statistics
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>% 
  group_by(gene, tissue) %>% 
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM)
