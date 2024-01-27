library(readxl)
library(dplyr)
library(tidyr)

excel_file_path <- 'Technical Test - Data Wrangling.xlsx'

# Specify sheet names
sheet_names <- c("RNA-seq (RPKM)", "Tissue Sample Metadata", "Serum Protein data", "Patient_clinical_data")

# Reading the given Excel file into a list of DataFrames
dfs_from_excel <- lapply(sheet_names, function(sheet) readxl::read_excel(excel_file_path, sheet = sheet, range = NULL, col_names = TRUE))

rnaseq <- dfs_from_excel[[1]]
tissue <- dfs_from_excel[[2]]
prot_meta <- dfs_from_excel[[3]]
patient <- dfs_from_excel[[4]]

# preprocessing
tissue <- tissue %>%
  rename(Patient_ID = `Patient  Number`, Sample_ID = Sample)
tissue$Patient_ID <- as.integer(tissue$Patient_ID)

prot_meta <- prot_meta %>%
  rename(Patient_ID = Patient, Sample_ID = Sample)
prot_meta$Patient_ID <- as.integer(prot_meta$Patient_ID)

patient <- patient %>%
  rename(Patient_ID = `Patient  Number`, Age = Age)
patient$Patient_ID <- as.integer(patient$Patient_ID)
patient$Age <- as.integer(patient$Age)

# Merge 'patient' and 'prot_meta' DataFrames
df <- patient %>%
  left_join(prot_meta, by = "Patient_ID", suffix = c("_left", "_right"))
df$Material_type <- 'SERUM'

# Melt the DataFrame to create a 'Gene_Symbol' column
merged_df1 <- df %>%
  gather(key = "Gene_Symbol", value = "Result", `Serum IL-6 (g/L)`, `Serum IL-6 Receptor (mg/L)`) %>%
  select(Study_ID, Patient_ID, Sex, Age, Sample_ID, Material_type, Gene_Symbol, Result)

# Map the gene symbols
gene_symbol_mapping <- c(`Serum IL-6 (g/L)` = 'IL6', `Serum IL-6 Receptor (mg/L)` = 'IL6R')
merged_df1$Gene_Symbol <- gene_symbol_mapping[merged_df1$Gene_Symbol]

merged_df1$Result_Units <- "g/L"
merged_df1$Sample_General_Pathology <- "NA"

merged_df1$Result <- as.numeric(as.character(merged_df1$Result))

# In the cases where we do not have a valid result, use "Not Done"
merged_df1$Status <- ifelse(is.na(merged_df1$Result) | !grepl("^\\d*\\.?\\d+$", as.character(merged_df1$Result)), 'NOT DONE', 'NA')

# take care of units
merged_df1 <- merged_df1 %>%
  mutate(Result = ifelse(Gene_Symbol == 'IL6R' & !is.na(Result), Result / 1000, Result),
         Status = ifelse(is.na(Result) | !grepl("^\\d*\\.?\\d+$", as.character(Result)), 'NOT DONE', 'NA'))


df2 <- patient %>%
  left_join(tissue, by = "Patient_ID") %>%
  select(-RIN, -`Total Reads(millions)`) %>%
  rename(Material_type = Material)

rnaseq_unpivoted <- rnaseq %>%
  gather(key = "Sample_ID", value = "Result", -GeneID)

# Merge 'df2' and 'rnaseq_unpivoted' DataFrames
merged_df2 <- df2 %>%
  left_join(rnaseq_unpivoted, by = "Sample_ID") %>%
  rename(Gene_Symbol = GeneID, Sample_General_Pathology = `Sample type`) %>%
  mutate(Status = ifelse(is.na(Result) | !grepl("^\\d*\\.?\\d+$", Result), 'NOT DONE', 'NA'),
         Result_Units = 'RPKM')

desired_columns <- c("Study_ID", "Patient_ID", "Sex", "Age", "Sample_ID", "Sample_General_Pathology", "Material_type", "Gene_Symbol", "Result", "Result_Units", "Status")

# Concatenate vertically 'merged_df2' and 'merged_df1'
example_report <- bind_rows(merged_df2[desired_columns], merged_df1[desired_columns])

example_report$Sex <- toupper(ifelse(example_report$Sex == 'M', 'Male', ifelse(example_report$Sex == 'F', 'Female', example_report$Sex)))

example_report$Unique_Patient_ID <- paste(example_report$Study_ID, example_report$Patient_ID, sep = '_')

# Reorder columns
example_report <- example_report %>%
  select(Study_ID, Patient_ID, Unique_Patient_ID, Sex, Age, Sample_ID, Sample_General_Pathology, Material_type, Gene_Symbol, Result, Result_Units, Status)

write.csv(example_report, file = 'example_report.csv', row.names = FALSE)
