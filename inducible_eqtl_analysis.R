#eQTL ANALYSIS PIPELINE####

#PACKAGES REQUIRED
library("edgeR")
library("csaw")
library("tidyverse")
library("dplyr")
library("tximport")
library("EnsDb.Hsapiens.v79")
library("gplots")
library("ggplot2")

#-----------------------------------------------------------------------------------------------------------------------------#

#01 EXCLUSION CRITERIA AND SAMPLE QUALITY CONTROL ####

#BRIEF: exclusion of patients meeting exclusion criteria

#file containing genotype, phenotype, ISARIC case report data for each patient
full_gt_phen_counts <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_phen_counts.csv")
#file containing ISARIC case report data
clinical_data <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/illness_severity.csv")

#EXCLUDE PATIENTS FOR WHOM GENOTYPE AT THE LEAD SNP LOCATION COULD NOT BE CALLED WITH GREATER THAN 90% PROBABILITY
full_gt_phen_counts <- full_gt_phen_counts %>% drop_na(Genotype)

#EXCLUDE UNRELATED INDIVIDUALS UP TO THE THIRD DEGREE
full_gt_phen_counts <- full_gt_phen_counts[!(full_gt_phen_counts$Unrelated == "False"),] #removes related individuals

#combine case report data
#get rid of duplicate values (data collected each day)
clinical_data <- clinical_data %>% dplyr::distinct(subjid, .keep_all = TRUE)
#create a vector of IDs for patients who remain included thus far and find case report data for these individuals 
IDs <- full_gt_phen_counts$canonical_isaric_id
clinical_data <- clinical_data[which(clinical_data$subjid %in% IDs), ]
#check step to see patients without case report data
str(setdiff(IDs,clinical_data$subjid))

#EXCLUDE PATIENTS WITH ISARIC CASE REPORT FORM UNAVAILABLE 
full_gt_phen_clinical <- merge(clinical_data, full_gt_phen_counts, by.x = 'subjid', by.y = 'canonical_isaric_id')

#save the file to ultra with a more accessible name to load in for downstream analysis
write.csv(full_gt_phen_clinical,"/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_phen_clinical.csv", row.names = F)
#-----------------------------------------------------------------------------------------------------------------------------#

#02 IMPORT SALMON TRANSCRIPT COUNTS####

#BRIEF: salmon transcript counts are imported

#file required for gene file names
full_gt_phen_clinical <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_phen_clinical.csv")

#obtain filenames of Salmon files
gene_count_file_names <- full_gt_phen_clinical$gene_count_file
file_names <- gsub(".salmon_gene_counts.csv$", "", gene_count_file_names)

#directory that contains Salmon files
rootpath <- "/scratch/u034/shared/wp5-rna-seq/analysis/rnaseq"
dirs <- dir(path = rootpath, recursive = TRUE, full.names = TRUE)

#to obtain Salmon files from remote directory
search <- dirs[grep(paste(file_names, collapse = "|"), dirs)]
salmon_files <- search[grep("^(?=.*quant.sf)(?!.*archive)", search, perl = TRUE)]

#check step to ensure files for all indivduals are available
file_names2 <- gsub(".quant.sf$", "", salmon_files)
setdiff(file_names2, file_names)

#to sort the vector in the same order as the id_gt_genecount
file_names <- file_names2[order(match(file_names2, file_names))]

#add quant.sf file at end of directory
salmon_files <- paste(file_names, "/quant.sf", sep = "")

#-----------------------------------------------------------------------------------------------------------------------------#

#03 SUMMARISE TRANSCRIPT COUNTS TO GENE LEVEL AND NORMALISE####

#BRIEF: salmon transcript counts are summarised to gene level and RNA-seq data normalised 

#load the annotation table for GrCh38
tx2gene <- read.csv("/scratch/u034/shared/wp5-rna-seq/analysis/rnaseq/20200627_NS550/salmon/tx2gene.csv")

files <- salmon_files
names(files) <- paste0(full_gt_phen_clinical$IID)

#IMPORT TRANSCRIPTS TO GENE LEVEL
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

normMat <- txi$length
cts <- txi$counts

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

#USE TMM NORMALISATION TO CORRECT FOR LIBRARY SIZE AND COMPOSITION
# Computing effective library sizes from scaled counts, to account for
# composition biases between samples
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

#to get group (genotype) for DGEList
genotype <- full_gt_phen_clinical$Genotype
group <- genotype

#creating a DGEList object
y <- DGEList(counts = cts, group = genotype)
y <- scaleOffset(y, normMat) 

#to filter out lowly expressed genes (remove untranscribed genes - those with 0 read counts)
keep <- filterByExpr(y)
y <- y[keep, ]

#to create a matrix of CPM
se <- SummarizedExperiment(assays = list(counts = y$counts, offset = y$offset))
se$totals <- y$samples$lib.size
cpms <- calculateCPM(se, use.offsets = TRUE, log = FALSE)
#LOG-TRANSFORM DATA
logcpm <- calculateCPM(se, use.offsets = TRUE, log = TRUE) 

#save the file to ultra with a more accessible name to load in for downstream analysis
logCPM_df <- as.data.frame(t(logcpm))
write.csv(logCPM_df,"/home/u034/jmulhollan/3p21.31_expression_analysis/logCPM_df.csv", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------#

#04 COMBINE GENOTYPE AND CLINICAL CHARACTERISTIC DATA WITH GENE COUNT DATA ####

#BRIEF: data containing genotype and clinical characteristic data is merged with the log normalised gene count data

full_gt_phen_clinical <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_phen_clinical.csv")
logCPM_df <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/logCPM_df.csv")

#merge the two datasets
full_gt_logCPM_phen_clinical <- cbind(full_gt_phen_clinical, logCPM_df)
str(full_gt_logCPM_phen_clinical)

#save the file to ultra with a more accessible name to load in for downstream analysis
write.csv(full_gt_logCPM_phen_clinical,"/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_logCPM_phen_clinical.csv", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------#

#05 LINEAR REGRESSION ####

#BRIEF: linear regression of gene expression on patient genotype 

#load in data
full_gt_logCPM_phen_clinical <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_logCPM_phen_clinical.csv")

#COVARIATES
#AGE
#data contains some 0 values and some NA values, 0 values may refer to infants
#test if an individual is an infant (<1 years old), if yes they should not have age changed
table(full_gt_logCPM_phen_clinical[full_gt_logCPM_phen_clinical$calc_age == 0, "subjid"], 
      full_gt_logCPM_phen_clinical[full_gt_logCPM_phen_clinical$calc_age == 0, "apdm_age"])
#replace NA values with mean age of full sample
full_gt_logCPM_phen_clinical$calc_age[is.na(full_gt_logCPM_phen_clinical$calc_age)] <- 
  mean(full_gt_logCPM_phen_clinical$calc_age, na.rm=TRUE)
#round mean values
full_gt_logCPM_phen_clinical$calc_age <- round(full_gt_logCPM_phen_clinical$calc_age, 0)

covariates <- c("Genotype", "Sex", "calc_age", "severity2", "sequencer_type", "PC1", "PC2", "PC3", "PC4", "PC5", 
                "PC6", "PC7", "PC8", "PC9", "PC10") #to include different covariates insert name here

#genes tested
genenames <- c("CCR2", "FYCO1", "CXCR6", "CCR3", 
               "LIMD1", "SACM1L", "SLC6A20",
               "LZTFL1", "CCR9", "XCR1", "CCR1") #to test different genes insert gene name here

#to convert to ensembl IDs
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genenames, keytype = "SYMBOL", columns = c("SYMBOL","GENEID")) 
genelist <- geneIDs$GENEID


#PERFORM LINEAR REGRESSION AND RETURN NORMALISED EFFECT SIZE AND P-VALUE FOR EACH GENE TESTED
glm_logCPM <- function(gene, ...) {
  combinedcovariates <- c(covariates, c(...)) 
  coef(summary(glm(as.formula(paste(gene,paste(combinedcovariates, collapse = " + "),
                                    sep = "~")), data = full_gt_logCPM_phen_clinical)))["Genotype",]
}

for (row in 1:nrow(geneIDs)) {
  Symbol <- geneIDs[row, "SYMBOL"]
  print(Symbol)
  Geneid <- geneIDs[row, "GENEID"]
  print(glm_logCPM(Geneid))
}

#-----------------------------------------------------------------------------------------------------------------------------#

#06 CALCULATE EXPRESSION CORRELATION AND PRODUCE HEATMAP ####

#BRIEF: calculates correlation for each gene relative to all other genes in the sample and plots it on a heatmap
#load in data
full_gt_logCPM_phen_clinical <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_logCPM_phen_clinical.csv")

genenames <- c("CCR2", "FYCO1", "CXCR6", "CCR3", 
               "LIMD1", "SACM1L", "SLC6A20",
               "LZTFL1", "CCR9", "XCR1", "CCR1") #to find correlation of different genes insert gene name here

#to convert to ensembl IDs
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genenames, keytype = "SYMBOL", columns = c("SYMBOL","GENEID")) 
genelist <- geneIDs$GENEID

#Create an empty vector to fill with the genes of interest IDs
genes_keep <- vector("character", nrow(geneIDs))
#loop to extract gene IDs
for (i in 1:nrow(geneIDs)) {
  x <- geneIDs$GENEID[i]
  genes_keep[i] <- x
}

#creates a dataframe of only the genes being investigated
keeps <- c(genes_keep) #isolates data from the stated columns
logCPM_counts <- full_gt_logCPM_phen_clinical[keeps] #dataframe of the isolated columns 

#to combine data and check correlation
logCPM_counts <- t(logCPM_counts)
logCPM_counts_bind<- cbind(geneIDs, logCPM_counts)
logCPM_counts_bind <- logCPM_counts_bind %>% remove_rownames %>% column_to_rownames(var="SYMBOL")
logCPM_counts_bind <- subset(logCPM_counts_bind, select = -c(GENEID))

#COMPUTING PEARSON'S CORRELATION COEFFICIENT
res <- cor(t(logCPM_counts_bind))

#PRODUCE A HEATMAP
col <- colorRampPalette(c("blue", "white", "red"))(20)
#formats order based on correlation
par(mar=c(5.1,4.1,4.1,2.1))
par(cex.main=1.5) #decreases title size

heatmap.2(x = res, col = col, margins=c(7,7), 
          density.info="none", trace="none", dendrogram='none',
          main = "Expression Correlation (log2CPM)",
          key.xlab = "Correlation of log2CPM", key.title = "Key", keysize = 1.4,
          cexRow = 1.4, cexCol = 1.4)

rpng.off()

#-----------------------------------------------------------------------------------------------------------------------------#

#07 GENE EXPRESSION FOR INDIVDUAL GENE ####

#BRIEF: view gene expression summary statistics for individual genes and plot a boxplot of expression

#load in data
full_gt_logCPM_phen_clinical <- read.csv("/home/u034/jmulhollan/3p21.31_expression_analysis/full_gt_logCPM_phen_clinical.csv")

#Enter gene in "______"
boxplotgene <- "CCR3" #manually input gene name here

#to find gene ID
geneID_name <- ensembldb::select(EnsDb.Hsapiens.v79, keys= boxplotgene, keytype = "SYMBOL", columns = c("SYMBOL","GENEID")) 
geneid <- geneID_name$GENEID

#create a dataframe of genotype and log2cpm counts for candidate genes
keeps <- c("Genotype", geneid) #isolates data from the stated columns
gene_expression_data <- full_gt_logCPM_phen_clinical[keeps] #dataframe of the isolated columns 

#need to change the names of Genotype to as they should appear on the boxplot
gene_expression_data$Genotype[gene_expression_data$Genotype == "0"] <- "Homozygous Major"
gene_expression_data$Genotype[gene_expression_data$Genotype == "1"] <- "Heterozygous"
gene_expression_data$Genotype[gene_expression_data$Genotype == "2"] <- "Homozygous Minor"
#need to make Genotype a factor 
gene_expression_data$Genotype <- factor(gene_expression_data$Genotype, levels = c("Homozygous Major", "Heterozygous", "Homozygous Minor"))

#mean and standard deviations 
homo_major <- gene_expression_data[which(gene_expression_data$Genotype=="Homozygous Major"),]
hetero <- gene_expression_data[which(gene_expression_data$Genotype=="Heterozygous"),]
homo_minor <- gene_expression_data[which(gene_expression_data$Genotype=="Homozygous Minor"),]

#homozygous major
mean(homo_major[,geneid])
sd(homo_major[,geneid])
#heterozygous
mean(hetero[,geneid])
sd(hetero[,geneid])
#homozygous minor
mean(homo_minor[,geneid])
sd(homo_minor[,geneid])

#creating a boxplot of the data
my_xlab <- paste(levels(gene_expression_data$Genotype),"\n(n=",table(gene_expression_data$Genotype),")",sep="") #no of observations

p <- ggplot(gene_expression_data, aes_string(x = "Genotype", y = geneid, fill="Genotype")) + 
  geom_boxplot() +
  labs(title = paste(boxplotgene, 'Gene Expression'), x = 'Genotype', y = paste(boxplotgene, 'normalised gene counts (log2CPM)')) + 
  stat_summary(fun=mean, geom="point", shape = 18, size = 4) 

p + theme(text = element_text(size=15), plot.title = element_text(hjust = 0.5, face = 'bold'), 
          axis.text=element_text(size=15), legend.position = "none") + 
  scale_fill_brewer(palette="Blues") + 
  scale_x_discrete(labels= my_xlab) #adds the number in each group below xlabs

#-----------------------------------------------------------------------------------------------------------------------------#
