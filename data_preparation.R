#DATA PREPARATION#

#PACKAGES####
library("tidyverse")
library("dplyr")
library("finalfit")


#FILES####

#BRIEF: Genotype data for ISARIC patients at the location of the SNP

#01 READING IN GENOTYPE DATA####
#phenotype data for ISARIC patients with genotyping data
raw_genotype_data <- read.table("geno.rs73064425.raw", header = TRUE)

#genotype column name is based on position (chr3.45859597.C.T_T) change this to genotype which is more accessible for downstream analyses
colnames(raw_genotype_data)[which(names(raw_genotype_data) == "chr3.45859597.C.T_T")] <- "Genotype"

#save the file to ultra with a more accessible name to load in for downstream analysis 
write.csv(raw_genotype_data, "/home/u034/bmedsci/3p21.31_expression_analysis/raw_genotype_data.csv", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------#

#02 READING IN PHENOTYPE DATA####

#BRIEF: Contains phenotype data used in analysis including ancestry

#phenotype data for ISARIC patientswith genotyping data
raw_phenotype_data <- read.csv("isaric.overlaprna.csv")

#save the file to ultra with a more accessible name to load in for downstream analysis
write.csv(raw_phenotype_data, "/home/u034/bmedsci/3p21.31_expression_analysis/raw_phenotype_data.csv", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------#

#03 READING IN RNA SEQUENCING DATA####

#BRIEF: Gene counts are imported directly from salmon, this file contains sequencer type

raw_gene_count_data <- read.csv("wp5_rnaseq_integration_20210208_120010.csv")

#visually inspect the sample ID for patients
raw_gene_count_data$isaric_sample_id 
#must be converted to lower case to match genotyping data
raw_gene_count_data$isaric_sample_id <- tolower(raw_gene_count_data$isaric_sample_id)

#save the file to ultra with a more accessible name to load in for downstream analysis
write.csv(raw_gene_count_data, "/home/u034/bmedsci/3p21.31_expression_analysis/raw_gene_count_data.csv", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------#

#04 MERGING GENOTYPE, PHENOTYPE, AND SEQUENCING DATA####

#BRIEF: combine genotype, phenotype and RNA seq data

#genotype
raw_genotype_data <- read.csv("/home/u034/bmedsci/3p21.31_expression_analysis/raw_genotype_data.csv")
#phenotype
raw_phenotype_data <- read.csv("/home/u034/bmedsci/3p21.31_expression_analysis/raw_phenotype_data.csv")
#RNA_seq
raw_gene_count_data <- read.csv("/home/u034/bmedsci/3p21.31_expression_analysis/raw_gene_count_data.csv")

#merge these dataframes centralising them to Isaric ID (IID)
gt_phen <- merge(raw_genotype_data, raw_phenotype_data, by.x = "IID", by.y = "X")
gt_phen_counts <- merge(gt_phen, raw_gene_count_data, by.x = "IID", by.y = "isaric_sample_id")
str(gt_phen_counts) #lets check the structure

#DEALING WITH DUPLICATES
str(gt_phen_counts) #there is 3 extra observations in  the merged dataset
setdiff(gt_phen$IID, gt_phen_counts$IID) #shows differences in values between two, as this is null must be duplicates
gt_phen_counts[duplicated(gt_phen_counts$IID),] #shows duplicated values in gen_phen_counts
#REMOVING DUPLICATED VALUES
full_gt_phen_counts <- gt_phen_counts[!duplicated(gt_phen_counts$IID),]

#save the file to ultra with a more accessible name to load in for downstream analysis
write.csv(full_gt_phen_counts, "/home/u034/bmedsci/3p21.31_expression_analysis/full_gt_phen_counts.csv", row.names = FALSE)

#to create a succinct dataframe containing ID, genotype and genecount files
keeps <- c("IID", "Genotype", "gene_count_file") #isolates data from the stated columns
ID_gt_genecounts <- full_gt_phen_counts[keeps] #dataframe of the isolated columns 
write.csv(ID_gt_genecounts, "/home/u034/bmedsci/3p21.31_expression_analysis/ID_gt_genecounts.csv", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------#

#05 ISARIC CASE REPORT DATA####

#BRIEF: create progression to ventilatory support or death variable in ISARIC case report data

#load in case report files, data is split between three files
ccp_data <- read.csv("/home/u034/bmedsci/3p21.31_expression_analysis/ccp_data.csv")
topline <- read.csv("/home/u034/bmedsci/3p21.31_expression_analysis/topline.csv")
outcome <- read.csv("/home/u034/bmedsci/3p21.31_expression_analysis/outcome.csv")

#join the death variable from topline
ccp_data = ccp_data %>%
  left_join(
    topline %>% dplyr::select(subjid, death), 
    by = c("subjid"))


# join the any_noninvasive, any_oxygen variable from outcome
ccp_data = ccp_data %>% 
  left_join(outcome %>% dplyr::select(subjid, any_noninvasive, any_oxygen), 
            by = c("subjid"))


#add the severity scale
ccp_data = ccp_data %>% 
  dplyr::mutate(
    severity2 = case_when(
      death == "Yes" ~ "Death",
      any_invasive == "Yes" ~ "IMV",
      any_noninvasive == "Yes" ~ "NIV/HFNC",
      daily_nasaloxy_cmtrt == "Yes" ~ "NIV/HFNC",
      any_oxygen == "Yes" ~ "Oxygen alone",
      TRUE ~ "Ward"
    ) %>% 
      factor(levels = c("Death", "IMV", "NIV/HFNC", "Oxygen alone", "Ward"))
  )

#this is a check step
ccp_data %>% 
  dplyr::select(subjid, death, any_invasive, any_noninvasive, daily_nasaloxy_cmtrt, any_oxygen, severity2)  %>% 
  head (20)

### fill the missing cestdat values 
ccp_data <- ccp_data %>% 
  group_by(subjid) %>%
  fill(cestdat) %>% 
  ungroup()

ccp_data <- as.data.frame(ccp_data)

#save the file to ultra with a more accessible name to load in for downstream analysis
write.csv(ccp_data,"/home/u034/bmedsci/3p21.31_expression_analysis/illness_severity.csv", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------#