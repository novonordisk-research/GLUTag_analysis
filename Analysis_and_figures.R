#### Load libraries ####

# General
library(tidyverse)

# Annotation
library(biomaRt)

# RNAseq data processing
library(tximeta)
library(sva)
library(DESeq2)

# Plotting
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplotify)
library(cowplot)
library(lemon)
library(ggpubr)

# For reporting
sessionInfo()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#### Gene conversion table ####

# Get a list of 1:1 homologues between mouse and human and respectives symbol IDs.

# gene ID conversion between species
#! Gria3 doesn't have orthologs between human and mouse in latest releases
mouse_biomart = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host ="https://nov2020.archive.ensembl.org") # latest mm10 is release 102 , cf. listEnsemblArchives() and http://www.ensembl.org/Mus_musculus/Info/Index
human_biomart = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host ="https://nov2020.archive.ensembl.org") # need to be on the same host
genes_conversion_table = getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), 
                                mart = mouse_biomart, 
                                attributesL = c("hgnc_symbol","ensembl_gene_id"), 
                                martL = human_biomart, 
                                uniqueRows= TRUE)
colnames(genes_conversion_table) <- c("mmusculus_ensembl_gene_id","mouse_symbol", "human_symbol","hsapiens_ensembl_gene_id")

# Restrict to 1:1 homologs
one2one_mouse <- genes_conversion_table %>% 
  group_by(mmusculus_ensembl_gene_id) %>% 
  summarise(n()) %>% 
  filter(`n()` <= 1) %>% 
  dplyr::select(mmusculus_ensembl_gene_id) %>% 
  pull()

one2one_human <- genes_conversion_table %>% 
  group_by(hsapiens_ensembl_gene_id) %>% 
  summarise(n()) %>% 
  filter(`n()` <= 1) %>% 
  dplyr::select(hsapiens_ensembl_gene_id) %>% 
  pull()

genes_conversion_table <- genes_conversion_table %>% 
  filter(mmusculus_ensembl_gene_id %in% one2one_mouse) %>%
  filter(hsapiens_ensembl_gene_id %in% one2one_human)

# Formatting
genes_conversion_table <- genes_conversion_table %>% 
  mutate(human_symbol = replace(human_symbol, human_symbol == "", NA))
genes_conversion_table <- genes_conversion_table %>% 
  mutate(mouse_symbol = replace(mouse_symbol, mouse_symbol == "", NA))

# Checks
head(genes_conversion_table)
dim(genes_conversion_table)
length(unique(genes_conversion_table$mmusculus_ensembl_gene_id))
length(unique(genes_conversion_table$hsapiens_ensembl_gene_id))
sum(duplicated(genes_conversion_table$mmusculus_ensembl_gene_id))
sum(duplicated(genes_conversion_table$hsapiens_ensembl_gene_id))
sum(duplicated(genes_conversion_table$mouse_symbol))
sum(duplicated(genes_conversion_table$human_symbol))
sum(is.na(genes_conversion_table$mouse_symbol))
sum(is.na(genes_conversion_table$human_symbol))


# Remove genes with no symbol IDs (which will be used for later data processing/filtering and analyses) or duplicated IDs

genes_conversion_table <- genes_conversion_table %>%
  filter(!is.na(mouse_symbol),
         !is.na(human_symbol),
         !duplicated(mouse_symbol),
         !duplicated(human_symbol))

# Checks
head(genes_conversion_table)
dim(genes_conversion_table)
length(unique(genes_conversion_table$mmusculus_ensembl_gene_id))
length(unique(genes_conversion_table$hsapiens_ensembl_gene_id))
sum(duplicated(genes_conversion_table$mmusculus_ensembl_gene_id))
sum(duplicated(genes_conversion_table$hsapiens_ensembl_gene_id))
sum(duplicated(genes_conversion_table$mouse_symbol))
sum(duplicated(genes_conversion_table$human_symbol))
sum(is.na(genes_conversion_table$mouse_symbol))
sum(is.na(genes_conversion_table$human_symbol))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

### CURRENTLY FIXING AN ISSUE IN THE FIRST SECTION BELOW (CHANGING THE CODE TO BE ABLE TO START FROM SRA FILE)

#### GSE193866: GLUTag mouse cell line  ####
# https://www.frontiersin.org/articles/10.3389/fphar.2022.861311/full

# Load sample metadata from former experiment
GLUTag_metadata <- read.csv(paste0("./GLUTag_new/SraRunTable.txt")) %>%
  select(Run,Sample.Name,Cell_Line,Cell_type,plate,PASSAGE)


# Extract MultiQC stats
GLUTag_metadata <- read.table("./GLUTag_new/Salmon_homemade/multiqc_data/multiqc_general_stats.txt",sep="\t",header = T) %>%
  as.data.frame() %>%
  filter(!grepl("Undetermined",Sample))  %>%
  dplyr::rename(Run = Sample) %>%
  mutate(Run = gsub(".*aux_info \\| ","",Run)) %>%
  rename_with(~gsub("Salmon_mqc.generalstats.salmon.","",.x)) %>%
  left_join(GLUTag_metadata,.) %>%
  dplyr::rename(Sample_Name = Sample.Name) %>%
  mutate(rownames = Run) %>%
  tibble::column_to_rownames("rownames") 

# Stats for paper
summary(GLUTag_metadata$percent_mapped)
summary(GLUTag_metadata$num_mapped)/1e6


# Load Salmon data
path_to_files <- "./GLUTag_new/Salmon_homemade"
files <- list.files(path=path_to_files, pattern="quant.sf", full.names = T, recursive = T)
names(files) <- files %>%
  gsub(path_to_files,"",.) %>%
  gsub("/|quant.sf","",.) 
GLUTag_gse <- tximeta(files,type = "salmon") %>%
  summarizeToGene(countsFromAbundance = "no")

# Generate count data
GLUTag_counts <- assays(GLUTag_gse)[["counts"]] %>%
  as.data.frame() %>%
  dplyr::select(matches(GLUTag_metadata$Run)) %>% 
  rownames_to_column("mmusculus_ensembl_gene_id") %>%
  mutate(mmusculus_ensembl_gene_id = gsub("\\..*","",mmusculus_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.) 
colnames(GLUTag_counts) <- gsub("_S.*_.*gz","",colnames(GLUTag_counts))

# Save count data for GEO
GLUTag_counts_GEO <- assays(GLUTag_gse)[["counts"]] %>%
  as.data.frame() %>%
  dplyr::select(matches(GLUTag_metadata$Run)) %>% 
  rownames_to_column("mmusculus_ensembl_gene_id")
colnames(GLUTag_counts_GEO) <- gsub("_S.*_.*gz","",colnames(GLUTag_counts_GEO))
head(GLUTag_counts_GEO)
write.table(GLUTag_counts_GEO,
            "./GLUTag_new/GEO_submission_Jan22/Raw_gene_counts_matrix.txt",
            row.names = F,
            col.names = T,
            quote = F)

# Checks
dim(GLUTag_counts)
head(GLUTag_counts)
sum(duplicated((GLUTag_counts$mmusculus_ensembl_gene_id)))

GLUTag_counts %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T) 
assays(GLUTag_gse)[["counts"]] %>% as.data.frame()  %>% dplyr::select(matches(GLUTag_metadata$Run)) %>% sum(na.rm=T)



# Generate CPM data
GLUTag_CPM <- assays(GLUTag_gse)[["counts"]] %>%
  as.data.frame() %>%
  rename_with(~gsub("_trimmed.fq.gz","",.x)) %>%
  dplyr::select(row.names(GLUTag_metadata)) %>% 
  as.data.frame() %>%
  mutate_all(funs((./sum(.))*1e6)) %>%
  as.data.frame() %>%
  dplyr::select(matches(GLUTag_metadata$Run)) %>% 
  rownames_to_column("mmusculus_ensembl_gene_id") %>%
  mutate(mmusculus_ensembl_gene_id = gsub("\\..*","",mmusculus_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.) 

# Checks
dim(GLUTag_CPM)
head(GLUTag_CPM)
sum(duplicated((GLUTag_CPM$mmusculus_ensembl_gene_id)))


# Checks
mean(colnames(GLUTag_counts[,-1][,-1]) == rownames(GLUTag_metadata))
mean(colnames(GLUTag_CPM[,-1][,-1]) == rownames(GLUTag_metadata))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#### GSE148224: Human intestinal organoid cell transcriptome ####
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148224

# Load sample metadata from the SraRunTable
GSE148224_metadata <- read.table("./GSE148224/SraRunTable.txt",sep=",",header=T)
GSE148224_metadata <- GSE148224_metadata %>%
  dplyr::select(Run,Sample.Name,Sample_number,facs_population,tissue,growth_medium)

# Extract MultiQC stats
GSE148224_metadata <- read.table("./GSE148224/Salmon_homemade_hg38_bootstraped/multiqc_data/multiqc_general_stats.txt",sep="\t",header = T) %>%
  as.data.frame() %>%
  filter(!grepl("Undetermined",Sample))  %>%
  dplyr::rename(Run = Sample) %>%
  mutate(Run = gsub(".*aux_info \\| ","",Run)) %>%
  rename_with(~gsub("Salmon_mqc.generalstats.salmon.","",.x)) %>%
  left_join(GSE148224_metadata,.) %>%
  dplyr::rename(Sample_Name = Sample.Name) %>%
  mutate(rownames = Run) %>%
  tibble::column_to_rownames("rownames")

# Stats for paper
summary(GSE148224_metadata$percent_mapped)
summary(GSE148224_metadata$num_mapped)/1e6


# Load Salmon data
path_to_files <- "./GSE148224/Salmon_homemade_hg38_bootstraped/"
files <- list.files(path=path_to_files, pattern="quant.sf", full.names = T, recursive = T)
names(files) <- files %>%
  gsub(path_to_files,"",.) %>%
  gsub("/|quant.sf","",.) 
GSE148224_gse <- tximeta(files,type = "salmon") %>%
  summarizeToGene(countsFromAbundance = "no")

# Generate count data
GSE148224_counts <- assays(GSE148224_gse)[["counts"]] %>%
  as.data.frame() %>%
  dplyr::select(matches(GSE148224_metadata$Run)) %>% 
  rownames_to_column("hsapiens_ensembl_gene_id") %>%
  mutate(hsapiens_ensembl_gene_id = gsub("\\..*","",hsapiens_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.) 

# Checks
dim(GSE148224_counts)
head(GSE148224_counts)
sum(duplicated((GSE148224_counts$mmusculus_ensembl_gene_id)))

GSE148224_counts %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)
assays(GSE148224_gse)[["counts"]] %>% as.data.frame()  %>% sum(na.rm=T)


# Generate TPM data
GSE148224_TPM <- assays(GSE148224_gse)[["abundance"]] %>%
  as.data.frame() %>%
  rename_with(~gsub("_trimmed.fq.gz","",.x)) %>%
  dplyr::select(matches(GSE148224_metadata$Run)) %>% 
  rownames_to_column("hsapiens_ensembl_gene_id") %>%
  mutate(hsapiens_ensembl_gene_id = gsub("\\..*","",hsapiens_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.)  


# checks
dim(GSE148224_TPM)
head(GSE148224_TPM)
sum(duplicated((GSE148224_TPM$mmusculus_ensembl_gene_id)))


# Checks
mean(colnames(GSE148224_counts[,-1][,-1]) == rownames(GSE148224_metadata))
mean(colnames(GSE148224_TPM[,-1][,-1]) == rownames(GSE148224_metadata))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#### GSE114913: Mouse duodenal enteroendocrine transcriptome ####
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114913

# Load sample metadata from the SraRunTable
GSE114913_metadata <- read.table("./GSE114913/SraRunTable.txt",sep=",",header=T)
GSE114913_metadata <- GSE114913_metadata %>%
  dplyr::select(Run,Sample.Name,cell_population,strain,tissue)
table(GSE114913_metadata$cell_population,GSE114913_metadata$strain)
table(GSE114913_metadata$Sample.Name,GSE114913_metadata$cell_population)
table(GSE114913_metadata$Sample.Name,GSE114913_metadata$strain)
table(GSE114913_metadata$Run,GSE114913_metadata$strain)

# Extract MultiQC stats
GSE114913_metadata <- read.table("./GSE114913/Salmon_homemade_mm10_bootstraped/multiqc_data/multiqc_general_stats.txt",sep="\t",header = T) %>%
  as.data.frame() %>%
  filter(!grepl("Undetermined",Sample))  %>%
  dplyr::rename(Run = Sample) %>%
  mutate(Run = gsub(".*aux_info \\| |_trim.*","",Run)) %>%
  rename_with(~gsub("Salmon_mqc.generalstats.salmon.","",.x)) %>%
  left_join(GSE114913_metadata,.) %>%
  dplyr::rename(Sample_Name = Sample.Name) %>%
  mutate(rownames = Run) %>%
  tibble::column_to_rownames("rownames")

# Stats
summary(GSE114913_metadata$percent_mapped)
summary(GSE114913_metadata$num_mapped)/1e6


# Load Salmon data
path_to_files <- "./GSE114913/Salmon_homemade_mm10_bootstraped/"
files <- list.files(path=path_to_files, pattern="quant.sf", full.names = T, recursive = T)
names(files) <- files %>%
  gsub(path_to_files,"",.) %>%
  gsub("/|quant.sf","",.) 
GSE114913_gse <- tximeta(files,type = "salmon") %>%
  summarizeToGene(countsFromAbundance = "no")


# Generate count data (and sum replicates together)
GSE114913_counts <- assays(GSE114913_gse)[["counts"]] %>%
  as.data.frame() %>%
  rename_with(~gsub("_trimmed.fq.gz","",.x)) %>%
  dplyr::select(matches(GSE114913_metadata$Run)) %>% 
  t()%>%
  as.data.frame()  %>%
  rownames_to_column("Run") %>%
  left_join(GSE114913_metadata[,c("Run","Sample_Name")],by="Run") %>%
  column_to_rownames("Run") %>%
  group_by(Sample_Name) %>%
  summarise_all(sum) %>% 
  column_to_rownames("Sample_Name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("mmusculus_ensembl_gene_id") %>%
  mutate(mmusculus_ensembl_gene_id = gsub("\\..*","",mmusculus_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.) 


# Checks
dim(GSE114913_counts)
head(GSE114913_counts)
sum(duplicated((GSE114913_counts$mmusculus_ensembl_gene_id)))

GSE114913_counts %>%  dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)
assays(GSE114913_gse)[["counts"]] %>% as.data.frame()  %>% sum(na.rm=T)

# Check before averaging TPM values for replicates
assays(GSE114913_gse)[["abundance"]] %>%
  as.data.frame() %>%
  rename_with(~gsub("_trimmed.fq.gz","",.x)) %>%
  dplyr::select(matches(GSE114913_metadata$Run)) %>%
  cor() %>%
  as.data.frame() %>%
  rownames_to_column("Run") %>%
  pivot_longer(-Run,names_to = "Run2", values_to = "TPM")  %>%
  left_join(GSE114913_metadata[,c("Run","Sample_Name")]) %>%
  dplyr::rename(Run1 = Run, 
                Sample_Name1 = Sample_Name,
                Run = Run2,
                Corr_TPM = TPM) %>%
  left_join(GSE114913_metadata[,c("Run","Sample_Name")]) %>%
  filter(Run1 != Run & Sample_Name1 == Sample_Name) %>%
  dplyr::select(Sample_Name,Corr_TPM) %>%
  distinct() %>%
  as.data.frame()
# Technical replicates of Neurod1-cre-eYFP samples are correlated > 0.9995: ok to average between replicates

# Generate TPM data (and average replicates together)
GSE114913_TPM <- assays(GSE114913_gse)[["abundance"]] %>%
  as.data.frame() %>%
  rename_with(~gsub("_trimmed.fq.gz","",.x)) %>%
  dplyr::select(matches(GSE114913_metadata$Run)) %>% 
  t()%>%
  as.data.frame()  %>%
  rownames_to_column("Run") %>%
  left_join(GSE114913_metadata[,c("Run","Sample_Name")],by="Run") %>%
  column_to_rownames("Run") %>%
  group_by(Sample_Name) %>%
  summarise_all(mean) %>% 
  column_to_rownames("Sample_Name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("mmusculus_ensembl_gene_id") %>%
  mutate(mmusculus_ensembl_gene_id = gsub("\\..*","",mmusculus_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.)  


# checks
dim(GSE114913_TPM)
head(GSE114913_TPM)
sum(duplicated((GSE114913_TPM$mmusculus_ensembl_gene_id)))



# Modify the metadata to reflect the merging of technical replicates
# Note: very similar percent_mapped and num_mapped values for the technical replicates
GSE114913_metadata <- GSE114913_metadata %>%
  dplyr::select(-Run) %>%
  group_by(Sample_Name) %>%
  summarise_at(vars(percent_mapped,num_mapped),mean) %>%
  left_join(GSE114913_metadata[,c("Sample_Name","cell_population","strain","tissue")]) %>%
  mutate(Run = Sample_Name) %>%
  dplyr::select(Run, everything(),percent_mapped,num_mapped) %>%
  distinct()
rownames(GSE114913_metadata) <- GSE114913_metadata$Run

# Checks
mean(colnames(GSE114913_counts[,-1][,-1]) == rownames(GSE114913_metadata))
mean(colnames(GSE114913_TPM[,-1][,-1]) == rownames(GSE114913_metadata))




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #





#### GSE114853: Human enteroendocrine cell transcriptome ####
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114853

# Load sample metadata from the SraRunTable
GSE114853_metadata <- read.table("./GSE114853/SraRunTable.txt",sep=",",header=T)
GSE114853_metadata <- GSE114853_metadata %>%
  dplyr::select(Run,Sample.Name,AvgSpotLen,source_name,facs_population,Participant)


# Extract MultiQC stats
GSE114853_metadata <- read.table("./GSE114853/Salmon_homemade_hg38_bootstraped/multiqc_data/multiqc_general_stats.txt",sep="\t",header = T) %>%
  as.data.frame() %>%
  filter(!grepl("Undetermined",Sample))  %>%
  dplyr::rename(Run = Sample) %>%
  mutate(Run = gsub(".*aux_info \\| |_trim.*","",Run)) %>%
  rename_with(~gsub("Salmon_mqc.generalstats.salmon.","",.x)) %>%
  left_join(GSE114853_metadata,.) %>%
  dplyr::rename(Sample_Name = Sample.Name) %>%
  mutate(rownames = Run) %>%
  tibble::column_to_rownames("rownames")

# Stats
summary(GSE114913_metadata$percent_mapped)
summary(GSE114913_metadata$num_mapped)/1e6


# Load Salmon data
path_to_files <- "./GSE114853/Salmon_homemade_hg38_bootstraped/"
files <- list.files(path=path_to_files, pattern="quant.sf", full.names = T, recursive = T)
names(files) <- files %>%
  gsub(path_to_files,"",.) %>%
  gsub("/|quant.sf","",.) 
GSE114853_gse <- tximeta(files,type = "salmon") %>%
  summarizeToGene(countsFromAbundance = "no")

# Generate count data
GSE114853_counts <- assays(GSE114853_gse)[["counts"]] %>%
  as.data.frame() %>%
  rename_with(~gsub("_trimmed.fq.gz","",.x)) %>%
  dplyr::select(matches(GSE114853_metadata$Run)) %>% 
  rownames_to_column("hsapiens_ensembl_gene_id") %>%
  mutate(hsapiens_ensembl_gene_id = gsub("\\..*","",hsapiens_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.) 


# Checks
dim(GSE114853_counts)
head(GSE114853_counts)
sum(duplicated((GSE114853_counts$mmusculus_ensembl_gene_id)))

GSE114853_counts %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)
assays(GSE114853_gse)[["counts"]] %>% as.data.frame() %>% sum(na.rm=T)


# Generate TPM data
GSE114853_TPM <- assays(GSE114853_gse)[["abundance"]] %>%
  as.data.frame() %>%
  rename_with(~gsub("_trimmed.fq.gz","",.x)) %>%
  dplyr::select(matches(GSE114853_metadata$Run)) %>% 
  rownames_to_column("hsapiens_ensembl_gene_id") %>%
  mutate(hsapiens_ensembl_gene_id = gsub("\\..*","",hsapiens_ensembl_gene_id)) %>%
  full_join(genes_conversion_table[,c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id")],.) 


# checks
dim(GSE114853_TPM)
head(GSE114853_TPM)
sum(duplicated((GSE114853_TPM$mmusculus_ensembl_gene_id)))


# Checks
mean(colnames(GSE114853_counts[,-1][,-1]) == rownames(GSE114853_metadata))
mean(colnames(GSE114853_TPM[,-1][,-1]) == rownames(GSE114853_metadata))





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#### Merge datasets ####

# Merge count data
merged_counts <- list(GLUTag_counts,GSE148224_counts,GSE114853_counts,GSE114913_counts) %>% 
  purrr::reduce(full_join, by = c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id"))
dim(merged_counts)

# Check that merging didn't introduce changes (duplicated rows etc) that would change the library size estimation
merged_counts %>% dplyr::select(matches(GLUTag_metadata$Sample_Name)) %>% sum(na.rm=T)
GLUTag_counts %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)

merged_counts %>% dplyr::select(matches(GSE148224_metadata$Run)) %>% sum(na.rm=T)
GSE148224_counts %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)

merged_counts %>% dplyr::select(matches(GSE114853_metadata$Run)) %>% sum(na.rm=T)
GSE114853_counts %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)

merged_counts %>% dplyr::select(matches(GSE114913_metadata$Run)) %>% sum(na.rm=T)
GSE114913_counts %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)


# Merge CPM/TPM data
merged_TPM <- list(GLUTag_CPM,GSE148224_TPM,GSE114853_TPM,GSE114913_TPM) %>% 
  purrr::reduce(full_join, by = c("hsapiens_ensembl_gene_id","mmusculus_ensembl_gene_id"))
dim(merged_TPM)

# Check that merging didn't introduce changes (duplicated rows etc) that would change the library size estimation
merged_TPM %>% dplyr::select(matches(GLUTag_metadata$Sample_Name)) %>% sum(na.rm=T)
GLUTag_CPM %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)

merged_TPM %>% dplyr::select(matches(GSE148224_metadata$Run)) %>% sum(na.rm=T)
GSE148224_TPM %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)

merged_TPM %>% dplyr::select(matches(GSE114853_metadata$Run)) %>% sum(na.rm=T)
GSE114853_TPM %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)

merged_TPM %>% dplyr::select(matches(GSE114913_metadata$Run)) %>% sum(na.rm=T)
GSE114913_TPM %>% dplyr::select(-hsapiens_ensembl_gene_id,-mmusculus_ensembl_gene_id) %>% sum(na.rm=T)



# Format and merge metadata
GLUTag_metadata$Dataset <- "GLUTag"
GSE148224_metadata$Dataset <- "GSE148224"
GSE114853_metadata$Dataset <- "GSE114853"
GSE114913_metadata$Dataset <- "GSE114913"
GLUTag_metadata$Dataset_2 <- "GLUTag"
GSE148224_metadata$Dataset_2 <- "Human_organoids"
GSE114853_metadata$Dataset_2 <- "Human_EECs"
GSE114913_metadata$Dataset_2 <- "Mouse_EECs"
merged_metadata <- do.call(plyr::rbind.fill,list(GLUTag_metadata,GSE148224_metadata,GSE114853_metadata,GSE114913_metadata))


merged_metadata <- merged_metadata %>%
  mutate(Condition = case_when(
    Dataset_2 == "GLUTag" ~ "GLUTag",
    Dataset_2 == "Human_EECs" ~ paste0(source_name,"_",facs_population),
    Dataset_2 == "Human_organoids" ~ facs_population,
    Dataset_2 == "Mouse_EECs" ~ paste0(strain,"_",cell_population) )) %>%
  mutate(Condition = gsub("GLP1-/CHGA-/SCG2-","non EEC",Condition)) %>%
  mutate(Condition = gsub("GLP1-/CHGA\\+/SCG2\\+","EEC GCG NEG",Condition)) %>%
  mutate(Condition = gsub("GLP1\\+/CHGA\\+/SCG2\\+","EEC GCG POS",Condition)) %>%
  mutate(Condition = gsub("GluVenus_negative","GLUVenus NEG",Condition)) %>%
  mutate(Condition = gsub("GluVenus_positive","GLUVenus POS",Condition)) %>%
  mutate(Condition = gsub("Neurod1-cre-eYFP_negative","NeuroD1 NEG",Condition)) %>%
  mutate(Condition = gsub("Neurod1-cre-eYFP_positive","NeuroD1 POS",Condition)) %>%
  mutate(Condition = gsub("Venus-","GLUVenus NEG",Condition)) %>%
  mutate(Condition = gsub("Venus\\+","GLUVenus POS",Condition)) %>%
  mutate(Condition = gsub(".*_","",Condition)) %>%
  mutate(Celltype = case_when(
    grepl("GLUVenus NEG",Condition) ~ "Mix of EEC and non-EEC, but no L-cells",
    grepl("non EEC|NeuroD1 NEG",Condition) ~ "Non-EEC cells",
    grepl("NeuroD1 POS",Condition) ~ "Total EECs population",
    grepl("EEC GCG NEG",Condition) ~ "Non-L-cell EECs",
    grepl("GLUVenus POS|GLUVenus POS|EEC GCG POS",Condition) ~ "L-cells",
    grepl("GLUTag",Condition) ~ Condition)) %>%
  mutate(Tissue = case_when(
    Dataset_2 == "Human_EECs" ~ source_name,
    Dataset_2 == "Human_organoids" ~ tissue ,
    Dataset_2 == "Mouse_EECs" ~ tissue,
    Dataset_2 == "GLUTag" ~ "cell line")) %>%
  mutate(Tissue = gsub("Ileal","ileal",Tissue)) %>%
  mutate(Tissue = gsub("Ileum","ileum",Tissue)) %>%
  mutate(Tissue = gsub("Jejunum","jejunum",Tissue)) %>%
  mutate(Species = case_when(
    Dataset_2 == "Human_EECs" ~ "Human",
    Dataset_2 == "Human_organoids" ~ "Human" ,
    Dataset_2 == "Mouse_EECs" ~ "Mouse",
    Dataset_2 == "GLUTag" ~ "Mouse")) %>%
  mutate(Condition_full = paste0(Species," ",Tissue,"\n",Condition)) %>%
  mutate(Tissue2 = paste(Species,Tissue))

table(merged_metadata$Celltype)
table(merged_metadata$Condition_full)
table(merged_metadata$Condition,merged_metadata$Celltype)
table(merged_metadata$Tissue,merged_metadata$Species)
table(merged_metadata$Tissue2)
table(merged_metadata$Condition_full ,merged_metadata$Dataset)

# Stats for paper
merged_metadata %>%
  group_by(Dataset) %>%
  summarise(mean_percent_mapped = mean(percent_mapped),
            sd_percent_mapped = sd(percent_mapped),
            mean_num_mapped = mean(num_mapped),
            sd_num_mapped = sd(num_mapped)) %>%
  mutate(mean_percent_mapped = round(mean_percent_mapped),
         sd_percent_mapped = round(sd_percent_mapped),
         mean_num_mapped = round(mean_num_mapped/1e6, digits=1),
         sd_num_mapped = round(sd_num_mapped/1e6, digits=1)
  )





# Checks
merged_counts <- merged_counts[,c("mmusculus_ensembl_gene_id","hsapiens_ensembl_gene_id",merged_metadata$Run)]
mean(colnames(merged_counts)[-1][-1] == merged_metadata$Run)

merged_TPM <- merged_TPM[,c("mmusculus_ensembl_gene_id","hsapiens_ensembl_gene_id",merged_metadata$Run)]
mean(colnames(merged_TPM)[-1][-1] == merged_metadata$Run)


mean(round(as.vector(colSums(merged_counts[-1][-1],na.rm=T))) == merged_metadata$num_mapped)
# differs for the datasets that have 2 technical replicates (ColSums will show 2X the number, as only the average is shown in num_mapped)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#### Batch correction, normalisation and filtering of count data ####

# Batch correction with ComBat_seq
merged_counts_adj <- ComBat_seq(as.matrix(merged_counts[,-1][,-1]), batch=factor(merged_metadata$Dataset))
merged_counts <- cbind(merged_counts[,c(1,2),drop=F],merged_counts_adj) %>%
  as.data.frame()


# VST transformation
# vst doesn't work with NA values. fix: replaced them with 0, then remove later on genes that add NAs
# (does this to avoid remove genes before normalisation by library size)
genes_to_include_post_vst <- merged_counts %>%
  drop_na() %>%
  pull(mmusculus_ensembl_gene_id)
merged_counts[is.na(merged_counts)] <- 0 

merged_counts <- cbind(DESeq2::vst(as.matrix(round(merged_counts[,-1][,-1]))),merged_counts[,c(1,2),drop=F])


# Gene filtering
# Remove genes that had NAs before vst()
# Remove genes with no symbol IDs
# Restrict to protein coding genes

gene_biotype <- getBM(attributes=c("mgi_symbol","gene_biotype"), mart=mouse_biomart)
protein_coding_list <- gene_biotype %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(mgi_symbol)

dim(merged_counts)
merged_counts <- merged_counts %>%
  left_join(genes_conversion_table,.) %>%
  filter(mmusculus_ensembl_gene_id %in% genes_to_include_post_vst) %>%
  dplyr::select(-mmusculus_ensembl_gene_id,-human_symbol,-hsapiens_ensembl_gene_id) %>%
  filter(!is.na(mouse_symbol)) %>%
  filter(mouse_symbol != "") %>%
  filter(mouse_symbol %in% protein_coding_list)
#filter(!grepl("mt-",mouse_symbol))
dim(merged_counts)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#### Extract gene lists ####
# From https://www.guidetopharmacology.org/download.jsp

targets_and_families <- read.csv("targets_and_families.csv")
table(targets_and_families$Type)


list_targets_and_families <- targets_and_families %>%
  dplyr::select(Family.name,MGI.symbol) %>%
  dplyr::rename(mouse_symbol = MGI.symbol) %>%
  filter(!is.na(mouse_symbol)) %>%
  filter(mouse_symbol != "") %>%
  pull(mouse_symbol)

list_gpcr_ion_channels <- targets_and_families %>%
  filter(Type %in% c("gpcr","vgic","lgic","other_ic")) %>%
  dplyr::select(Family.name,MGI.symbol) %>%
  dplyr::rename(mouse_symbol = MGI.symbol) %>%
  filter(!is.na(mouse_symbol)) %>%
  filter(mouse_symbol != "") %>%
  pull(mouse_symbol)


list_glutamatergic_receptors <- read.csv("targets_and_families.csv") %>%
  filter(grepl("glutamat",Family.name)) %>%
  dplyr::select(Family.name,MGI.symbol) %>%
  dplyr::rename(mouse_symbol = MGI.symbol) %>%
  filter(!is.na(mouse_symbol)) %>%
  filter(mouse_symbol != "") %>%
  filter(grepl("receptors",Family.name)) %>%
  mutate(Receptor_type = case_when(
    grepl("Grm",mouse_symbol) ~ "Metabotropic receptors",
    grepl("Gria",mouse_symbol) ~ "Ionotropic receptors (AMPA)",
    grepl("Grin",mouse_symbol) ~ "Ionotropic receptors (NMDA)",
    grepl("Grid",mouse_symbol) ~ "Ionotropic receptors",
    grepl("Grik",mouse_symbol) ~ "Ionotropic receptors (Kainate)"
  ))

table(list_glutamatergic_receptors$Receptor_type)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### Save processed data ####

save.image(file='Paper_code_Environment.RData')
#load('Paper_code_Environment.RData')




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### Annotations / color settings for figures ####

annotation_col <-  merged_metadata %>%
  dplyr::select(Condition,Celltype,Tissue,Species,Condition_full,Tissue2) %>%
  mutate_all(as.factor) %>%
  as_tibble() %>%
  dplyr::select(Condition_full,Celltype,Tissue2) %>%
  dplyr::rename(`Cell type` = Celltype) %>%
  distinct() %>%
  column_to_rownames("Condition_full") %>%
  dplyr::rename(Tissue = Tissue2) %>%
  mutate(`Cell type` = factor(`Cell type`, levels = c("GLUTag","L-cells","Non-L-cell EECs","Total EECs population","Non-EEC cells","Mix of EEC and non-EEC, but no L-cells")),
         Tissue = factor(Tissue, levels = c("Mouse cell line","Mouse duodenum","Human jejunum","Human ileum","Human ileal organoids")))

Celltype_colours <- brewer.pal(9,"Set1")[1:6]
names(Celltype_colours) <- levels(annotation_col$`Cell type`)
Tissue_colours <- brewer.pal(8,"Set2")[1:5]
names(Tissue_colours) <- levels(annotation_col$Tissue)
Tissue_colours["Mouse cell line"] <- Celltype_colours["GLUTag"]

ann_colors = list(
  `Cell type` = Celltype_colours,
  Tissue = Tissue_colours
)

Tissue_colours_geneplot <- Tissue_colours
names(Tissue_colours_geneplot) <- gsub("Mouse cell line","GLUTag",names(Tissue_colours_geneplot))
names(Tissue_colours_geneplot) <- gsub(" ","\n",names(Tissue_colours_geneplot))


#### Sample-level plots - PCA ####

# Perform PCA in the batch-corrected vst-transformed data
exp_raw <- merged_counts %>%
  column_to_rownames("mouse_symbol")
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

# Merge the PCA results with the metadata
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],merged_metadata[,c("Run","Tissue2","Celltype")]) %>%
  mutate(`Cell type` = factor(Celltype, levels = c("GLUTag","L-cells","Non-L-cell EECs","Total EECs population","Non-EEC cells","Mix of EEC and non-EEC, but no L-cells")),
         Tissue = factor(Tissue2, levels = c("Mouse cell line","Mouse duodenum","Human jejunum","Human ileum","Human ileal organoids")))

# Generate PCA plot
PCA_plot <- ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(fill = `Cell type`,colour = `Cell type`, shape = Tissue), size=4, colour="black") +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_classic()  +
  scale_fill_manual(breaks=names(Celltype_colours),values=Celltype_colours) +
  scale_colour_manual(breaks=names(Celltype_colours),values=Celltype_colours) +
  scale_shape_manual(breaks=names(Tissue_colours),values=c(21:25)) +
  theme(text = element_text(size=20)) +
  theme(strip.text.y.right = element_text(angle = 0))  +
  guides(fill=guide_legend(override.aes=list(shape=22))) + 
  guides(shape = guide_legend(order = 1),
         colour = guide_legend(order = 2))
PCA_plot


# Calculate the contribution of each genes
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#pca-results-for-variables
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
loadings <- PCA_raw$rotation
sdev <- PCA_raw$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))

# Extract the top 20 genes for the first PC
PC1_genes <- var.contrib %>%
  as.data.frame() %>%
  arrange(-abs(PC1)) %>%
  rownames_to_column("symbol") %>%
  dplyr::select(symbol,PC1) %>%
  mutate(PC1 = round(PC1,2)) %>%
  mutate(PC1 = paste0("(",PC1,"%)")) %>%
  dplyr::select(symbol,PC1) %>%
  as.data.frame() %>%
  head(20)

# Extract the top 20 genes for the second PC
PC2_genes <- var.contrib %>%
  as.data.frame() %>%
  arrange(-abs(PC2)) %>%
  rownames_to_column("symbol") %>%
  dplyr::select(symbol,PC2) %>%
  mutate(PC2 = round(PC2,2)) %>%
  mutate(PC2 = paste0("(",PC2,"%)")) %>%
  dplyr::select(symbol,PC2) %>%
  as.data.frame() %>%
  head(20)

# Create a result table
PC_genes <- cbind(seq(1:length(PC1_genes$symbol)), PC1_genes, PC2_genes)
colnames(PC_genes) <- c("#","Top genes: PC1","","Top genes: PC2","")
Table_features <- ggtexttable(PC_genes, rows = NULL, theme = ttheme("blank", base_size = 18)) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>%
  table_cell_font(column = c(3,5), face = "italic", row=2:21, size=16)
Table_features


#### Sample-level plots - heatmaps ####

# Heatmap for all genes

N_genes <- dim(merged_counts)[1]
cor_table <- merged_counts %>%
  pivot_longer(-mouse_symbol, names_to = "Run", values_to = "Counts") %>%
  left_join(merged_metadata[,c("Run","Condition_full")]) %>%
  group_by(mouse_symbol, Condition_full) %>%
  summarise(mean = mean(Counts,na.rm=T)) %>%
  pivot_wider(names_from = "Condition_full",values_from = "mean") %>%
  column_to_rownames("mouse_symbol") %>%
  cor(use = "complete.obs", method="pearson")

#check there is not negative correlation before squaring
mean(cor_table>0)
cor_table = cor_table ^ 2

N_cut = 2
annotation_col <- annotation_col[match(colnames(cor_table),row.names(annotation_col)),]
heatmap_all <- cor_table %>%
  as.data.frame() %>%
  rownames_to_column("rownames") %>%
  mutate(rownames = gsub("\n"," ",rownames)) %>%
  column_to_rownames("rownames") %>%
  ComplexHeatmap::pheatmap(., 
                           display_numbers = TRUE, 
                           annotation_col = annotation_col,
                           annotation_colors = ann_colors,
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean", 
                           clustering_method = "complete",
                           #color =  circlize::colorRamp2(c(min(.),1), c("white", "darkorange2")),
                           angle_col = "45",
                           treeheight_col = 20,
                           name = "Pearson R2",
                           scale = "none",
                           border_color = "black",
                           cutree_rows = N_cut,
                           cutree_cols = N_cut,
                           fontsize = 8,
                           main = paste0("Correlation heatmap (Protein-coding genes, n=",N_genes,")"))


# Dendrogram only
dendro_all <- ComplexHeatmap::pheatmap(cor_table, 
                                       color = c("white","white"),
                                       annotation_row = annotation_col,
                                       annotation_colors = ann_colors,
                                       angle_col = "45",
                                       name = " ",
                                       scale = "none",
                                       cellwidth = 0,
                                       treeheight_col = 0,
                                       labels_col = 0,
                                       show_colnames = F,
                                       border_color = NA,
                                       legend=F,
                                       annotation_legend = F,
                                       main =paste0( " Protein-coding genes (n=",N_genes,")\n\n\n")) 



# Most correlated samples table
cor_table2 <- cor_table
diag(cor_table2)<-NA
colnames(cor_table2) <- gsub("\\\n"," ",colnames(cor_table2))
row.names(cor_table2) <- gsub("\\\n"," ",row.names(cor_table2))
top_cors <- apply(cor_table2,1,function(x){paste0(rownames(cor_table2)[order(x, decreasing =TRUE)][1:4],
                                                  " (",
                                                  round(x[order(x, decreasing =TRUE)][1:4],digits=3),
                                                  ")")}) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("1st most similar", "2nd most similar", "3rd most similar", "4th most similar")) %>%
  rownames_to_column("Sample") 
top_cors





# Heatmap for the "target and families" genes
N_genes <- merged_counts %>%
  filter(mouse_symbol %in% list_targets_and_families) %>%
  pull(mouse_symbol) %>% 
  length()

cor_table <- merged_counts %>%
  filter(mouse_symbol %in% list_targets_and_families) %>%
  pivot_longer(-mouse_symbol, names_to = "Run", values_to = "Counts") %>%
  left_join(merged_metadata[,c("Run","Condition_full")]) %>%
  group_by(mouse_symbol, Condition_full) %>%
  summarise(mean = mean(Counts,na.rm=T)) %>%
  pivot_wider(names_from = "Condition_full",values_from = "mean") %>%
  column_to_rownames("mouse_symbol") %>%
  cor(use = "complete.obs", method="pearson") 

#check there is not negative correlation before squaring
mean(cor_table>0)
cor_table = cor_table ^ 2



N_cut = 2
annotation_col <- annotation_col[match(colnames(cor_table),row.names(annotation_col)),]
heatmap_ligandsreceptors <- cor_table %>%
  as.data.frame() %>%
  rownames_to_column("rownames") %>%
  mutate(rownames = gsub("\n"," ",rownames)) %>%
  column_to_rownames("rownames") %>%
  ComplexHeatmap::pheatmap(.,  
                           display_numbers = TRUE, 
                           annotation_col = annotation_col,
                           annotation_colors = ann_colors,
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean", 
                           clustering_method = "complete",
                           #color =  circlize::colorRamp2(c(min(.),1), c("white", "darkorange2")),
                           angle_col = "45",
                           treeheight_col = 20,
                           name = "Pearson R2",
                           scale = "none",
                           border_color = "black",
                           cutree_rows = N_cut,
                           cutree_cols = N_cut,
                           fontsize = 8,
                           main = paste0("Correlation heatmap (Ligands & target proteins genes, n=",N_genes,")"))


# Dendrogram only
dendro_ligandsreceptors <- ComplexHeatmap::pheatmap(cor_table, 
                                                    color = c("white","white"),
                                                    annotation_row = annotation_col,
                                                    annotation_colors = ann_colors,
                                                    angle_col = "45",
                                                    name = "  ",
                                                    scale = "none",
                                                    cellwidth = 0,
                                                    treeheight_col = 0,
                                                    labels_col = 0,
                                                    show_colnames = F,
                                                    border_color = NA,
                                                    #border_color = "black",
                                                    legend=F,
                                                    annotation_legend = F,
                                                    main = paste0(" Ligands & target proteins genes (n=",N_genes,")\n\n\n")) 




# Most correlated samples table
cor_table2 <- cor_table
diag(cor_table2)<-NA
colnames(cor_table2) <- gsub("\\\n"," ",colnames(cor_table2))
row.names(cor_table2) <- gsub("\\\n"," ",row.names(cor_table2))
top_cors <- apply(cor_table2,1,function(x){paste0(rownames(cor_table2)[order(x, decreasing =TRUE)][1:4],
                                                  " (",
                                                  round(x[order(x, decreasing =TRUE)][1:4],digits=3),
                                                  ")")}) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("1st most similar", "2nd most similar", "3rd most similar", "4th most similar")) %>%
  rownames_to_column("Sample") 
top_cors






# Heatmap for the "gpcr & ion channels" genes
N_genes <- merged_counts %>%
  filter(mouse_symbol %in% list_gpcr_ion_channels) %>%
  pull(mouse_symbol) %>% 
  length()

cor_table <- merged_counts %>%
  filter(mouse_symbol %in% list_gpcr_ion_channels) %>%
  pivot_longer(-mouse_symbol, names_to = "Run", values_to = "Counts") %>%
  left_join(merged_metadata[,c("Run","Condition_full")]) %>%
  group_by(mouse_symbol, Condition_full) %>%
  summarise(mean = mean(Counts,na.rm=T)) %>%
  pivot_wider(names_from = "Condition_full",values_from = "mean") %>%
  column_to_rownames("mouse_symbol") %>%
  cor(use = "complete.obs", method="pearson") 

#check there is not negative correlation before squaring
mean(cor_table>0)
cor_table = cor_table ^ 2




N_cut = 2
annotation_col <- annotation_col[match(colnames(cor_table),row.names(annotation_col)),]
heatmap_gpcrionchannel <- cor_table %>%
  as.data.frame() %>%
  rownames_to_column("rownames") %>%
  mutate(rownames = gsub("\n"," ",rownames)) %>%
  column_to_rownames("rownames") %>%
  ComplexHeatmap::pheatmap(., 
                           display_numbers = TRUE, 
                           annotation_col = annotation_col,
                           annotation_colors = ann_colors,
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean", 
                           clustering_method = "complete",
                           #color =  circlize::colorRamp2(c(min(.),1), c("white", "darkorange2")),
                           angle_col = "45",
                           treeheight_col = 20,
                           name = "Pearson R2",
                           scale = "none",
                           border_color = "black",
                           cutree_rows = N_cut,
                           cutree_cols = N_cut,
                           fontsize = 8,
                           main = paste0("Correlation heatmap (GPCRs & ion channels genes, n=",N_genes,")"))



# Dendrogram only
dendro_gpcrionchannel <- ComplexHeatmap::pheatmap(cor_table, 
                                                  color = c("white","white"),
                                                  annotation_row = annotation_col,
                                                  annotation_colors = ann_colors,
                                                  angle_col = "45",
                                                  name = "   ",
                                                  scale = "none",
                                                  cellwidth = 0,
                                                  treeheight_col = 0,
                                                  labels_col = 0,
                                                  show_colnames = F,
                                                  border_color = NA,
                                                  #cutree_rows = N_cut,
                                                  legend=F,
                                                  main = paste0(" GPCRs & ion channels genes (n=",N_genes,")\n\n\n"))






# Most correlated samples table
cor_table2 <- cor_table
diag(cor_table2)<-NA
colnames(cor_table2) <- gsub("\\\n"," ",colnames(cor_table2))
row.names(cor_table2) <- gsub("\\\n"," ",row.names(cor_table2))
top_cors <- apply(cor_table2,1,function(x){paste0(rownames(cor_table2)[order(x, decreasing =TRUE)][1:4],
                                                  " (",
                                                  round(x[order(x, decreasing =TRUE)][1:4],digits=3),
                                                  ")")}) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("1st most similar", "2nd most similar", "3rd most similar", "4th most similar")) %>%
  rownames_to_column("Sample") 
top_cors












# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#### Gene-level plots ####

Expression_plot_data <- merged_TPM %>%
  dplyr::select(-hsapiens_ensembl_gene_id) %>%
  as.data.frame() %>%
  dplyr::select(mmusculus_ensembl_gene_id, matches(merged_metadata$Run)) %>% 
  full_join(genes_conversion_table[,c("mouse_symbol","mmusculus_ensembl_gene_id")],.) %>%
  left_join(list_glutamatergic_receptors,.) %>%
  dplyr::select(-mmusculus_ensembl_gene_id,-Family.name) %>%
  pivot_longer(-c(Receptor_type, mouse_symbol), names_to = "Run", values_to = "TPM") %>%
  left_join(merged_metadata)
table(Expression_plot_data$Receptor_type)


list_geneplots <- list()

for (subset_name in c("NMDA","Metabo","AMPA","Kaina")){
  
  print(subset_name)
  
  
  if (subset_name == "NMDA"){
    
    list_geneplots[[paste0(subset_name,"_1")]] <- Expression_plot_data %>%
      filter(Celltype %in% c("GLUTag","L-cells")) %>%
      filter(grepl(subset_name,Receptor_type)) %>%
      mutate(Condition_Lcells = gsub("\n.*","",Condition_full)) %>%
      mutate(Condition_Lcells = gsub("Mouse cell line","GLUTag",Condition_Lcells)) %>%
      mutate(Condition_Lcells = gsub(" ","\n",Condition_Lcells)) %>%
      mutate(Condition_Lcells = factor(Condition_Lcells,levels=c("GLUTag","Mouse\nduodenum","Human\njejunum","Human\nileum","Human\nileal\norganoids"))) %>%
      ggplot(aes(x=mouse_symbol,y=log(TPM+1), fill=Condition_Lcells)) +
      geom_boxplot(outlier.shape = NA) + 
      guides(fill = "none") +
      theme_classic() +
      geom_jitter(size=1) +
      xlab(" ") +
      facet_rep_grid(Condition_Lcells~.,scales = "free", repeat.tick.labels = T) +
      scale_fill_manual(breaks=names(Tissue_colours_geneplot),values=Tissue_colours_geneplot) +
      #ggtitle("Expression levels\n") + 
      theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
            text = element_text(size=14)) +
      theme(strip.text.y.right = element_text(angle = 0))
    
    list_geneplots[[paste0(subset_name,"_2")]] <- Expression_plot_data %>%
      filter(Celltype %in% c("GLUTag","L-cells")) %>%
      filter(grepl(subset_name,Receptor_type)) %>%
      mutate(TPM = log(TPM+1)) %>%
      dplyr::select(mouse_symbol,Condition_full, TPM) %>%
      group_by(Condition_full,mouse_symbol) %>%
      summarise(TPM_mean = mean(TPM)) %>%
      mutate(Condition_Lcells = gsub("\n.*","",Condition_full)) %>%
      mutate(Condition_Lcells = gsub("Mouse cell line","GLUTag",Condition_Lcells)) %>%
      mutate(Condition_Lcells = gsub(" ","\n",Condition_Lcells)) %>%
      mutate(Condition_Lcells = factor(Condition_Lcells,levels=c("GLUTag","Mouse\nduodenum","Human\njejunum","Human\nileum","Human\nileal\norganoids"))) %>%
      ungroup() %>%
      dplyr::select(-Condition_full) %>%
      pivot_wider(names_from = "Condition_Lcells", values_from = "TPM_mean") %>%
      column_to_rownames("mouse_symbol") %>%
      #mutate_all(~scale(.,center=F, scale=T)) %>%
      mutate_all(~scales::rescale(.,to = c(0, 1))) %>%
      ComplexHeatmap::pheatmap(.,
                               clustering_distance_cols = "euclidean", 
                               clustering_method = "complete",
                               color =  circlize::colorRamp2(c(0,1), c("white", "red")),
                               scale = "none", 
                               name = "Rescaled\nexpression\nper sample",
                               drop_levels = TRUE,
                               border_color = "black",
                               cluster_rows = F,
                               # cellheight = 30,
                               # cellwidth = 60,
                               # angle_col = "0",
                               fontsize = 12)#,
    #main = "Relative expression levels\n")
    
    cols <- rep("black",dim(list_geneplots[[paste0(subset_name,"_2")]]@matrix)[2])
    cols[which(colnames(list_geneplots[[paste0(subset_name,"_2")]]@matrix) ==  "GLUTag")] <- "red"
    list_geneplots[[paste0(subset_name,"_2")]]@column_names_param$gp <- grid::gpar(col=cols, fontsize=12)
    
    
    list_geneplots[[paste0(subset_name,"_2")]] <- grid::grid.grabExpr(draw(list_geneplots[[paste0(subset_name,"_2")]]))
    
  } else {
    
    list_geneplots[[paste0(subset_name,"_1")]] <- Expression_plot_data %>%
      filter(Celltype %in% c("GLUTag","L-cells")) %>%
      filter(grepl(subset_name,Receptor_type)) %>%
      mutate(Condition_Lcells = gsub("\n.*","",Condition_full)) %>%
      mutate(Condition_Lcells = gsub("Mouse cell line","GLUTag",Condition_Lcells)) %>%
      mutate(Condition_Lcells = gsub(" ","\n",Condition_Lcells)) %>%
      mutate(Condition_Lcells = factor(Condition_Lcells,levels=c("GLUTag","Mouse\nduodenum","Human\njejunum","Human\nileum","Human\nileal\norganoids"))) %>%
      ggplot(aes(x=mouse_symbol,y=log(TPM+1), fill=Condition_Lcells)) +
      geom_boxplot(outlier.shape = NA) + 
      guides(fill = "none") +
      theme_classic() +
      geom_jitter(size=0.5) +
      xlab(" ") +
      facet_rep_grid(Condition_Lcells~.,scales = "free", repeat.tick.labels = T) +
      scale_fill_manual(breaks=names(Tissue_colours_geneplot),values=Tissue_colours_geneplot) +
      #ggtitle("Expression levels\n") + 
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            text = element_text(size=11)) +
      theme(strip.text.y.right = element_text(angle = 0))
    
    
    list_geneplots[[paste0(subset_name,"_2")]] <- Expression_plot_data %>%
      filter(Celltype %in% c("GLUTag","L-cells")) %>%
      filter(grepl(subset_name,Receptor_type)) %>%
      mutate(TPM = log(TPM+1)) %>%
      dplyr::select(mouse_symbol,Condition_full, TPM) %>%
      group_by(Condition_full,mouse_symbol) %>%
      summarise(TPM_mean = mean(TPM)) %>%
      mutate(Condition_Lcells = gsub("\n.*","",Condition_full)) %>%
      mutate(Condition_Lcells = gsub("Mouse cell line","GLUTag",Condition_Lcells)) %>%
      mutate(Condition_Lcells = gsub(" ","\n",Condition_Lcells)) %>%
      mutate(Condition_Lcells = factor(Condition_Lcells,levels=c("GLUTag","Mouse\nduodenum","Human\njejunum","Human\nileum","Human\nileal\norganoids"))) %>%
      ungroup() %>%
      dplyr::select(-Condition_full) %>%
      pivot_wider(names_from = "Condition_Lcells", values_from = "TPM_mean") %>%
      column_to_rownames("mouse_symbol") %>%
      #mutate_all(~scale(.,center=F, scale=T)) %>%
      mutate_all(~scales::rescale(.,to = c(0, 1))) %>%
      ComplexHeatmap::pheatmap(.,
                               clustering_distance_cols = "euclidean", 
                               clustering_method = "complete",
                               color =  circlize::colorRamp2(c(0,1), c("white", "red")),
                               scale = "none", 
                               name = "Rescaled\nexpression\nper sample",
                               drop_levels = TRUE,
                               border_color = "black",
                               cluster_rows = F,
                               # cellheight = 30,
                               # cellwidth = 60,
                               # angle_col = "0",
                               fontsize = 10)#,
    #main = "Relative levels\n")
    
    
    # write the GLUTag rownames in red
    
    cols <- rep("black",dim(list_geneplots[[paste0(subset_name,"_2")]]@matrix)[2])
    cols[which(colnames(list_geneplots[[paste0(subset_name,"_2")]]@matrix) ==  "GLUTag")] <- "red"
    list_geneplots[[paste0(subset_name,"_2")]]@column_names_param$gp <- grid::gpar(col=cols, fontsize=10)
    
    list_geneplots[[paste0(subset_name,"_2")]] <- grid::grid.grabExpr(draw(list_geneplots[[paste0(subset_name,"_2")]]))
    
  }
  
}




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#### Figures for the paper ####

###### Figure 1 ######

# write the GLUTag rownames in red
cols <- rep("black",dim(dendro_all@matrix)[1])
cols[which(row.names(dendro_all@matrix) ==  "Mouse cell line\nGLUTag")] <- "red"
dendro_all@row_names_param$gp <- grid::gpar(col=cols, fontsize=10)

cols <- rep("black",dim(dendro_ligandsreceptors@matrix)[1])
cols[which(row.names(dendro_ligandsreceptors@matrix) ==  "Mouse cell line\nGLUTag")] <- "red"
dendro_ligandsreceptors@row_names_param$gp <- grid::gpar(col=cols, fontsize=10)

cols <- rep("black",dim(dendro_gpcrionchannel@matrix)[1])
cols[which(row.names(dendro_gpcrionchannel@matrix) ==  "Mouse cell line\nGLUTag")] <- "red"
dendro_gpcrionchannel@row_names_param$gp <- grid::gpar(col=cols, fontsize=10)


# Convert heatmap objects to a format compatible with cowplot
dendro_all <- grid::grid.grabExpr(draw(dendro_all))
dendro_ligandsreceptors <- grid::grid.grabExpr(draw(dendro_ligandsreceptors))
dendro_gpcrionchannel <- grid::grid.grabExpr(draw(dendro_gpcrionchannel))


# Set titles
title_top <- ggdraw() + 
  draw_label(
    "Sample clustering (complete linkage hierarchical clustering using the Euclidean distances)",
    fontface = 'bold',
    hjust = 0.5,
    size = 16
  )

title_bottom <- ggdraw() + 
  draw_label(
    "NMDA receptor subunits expression in L-cells and related cell models",
    fontface = 'bold',
    hjust = 0.5,
    size = 16
  )





# Generate plot 
top_row = plot_grid(dendro_all,dendro_ligandsreceptors,NULL,dendro_gpcrionchannel,rel_widths = c(0.7,0.7,0.1,1), nrow=1)
bottom_row = plot_grid(list_geneplots[["NMDA_1"]],NULL,list_geneplots[["NMDA_2"]],rel_widths = c(1, 0.05, 1), nrow=1)
plot_grid(title_top,
          top_row,
          NULL,
          title_bottom,
          bottom_row,
          rel_heights = c(0.1, 0.9, 0.05, 0.1, 0.9), 
          nrow = 5,
          labels = c("A"," "," ","B"," "," "))

# Only add the lines needed (these are not set beforehand, else black lines would appears for the heatmap even if tits dimension has been set to 0)
list_to_edit <- grid.ls(grid.force(),flatten = T, print=F) %>%
  unlist() %>%
  as.data.frame() %>%
  filter(grepl("GRID.rect",.)) %>%
  pull(.)
list_to_edit <- list_to_edit[-c(1,2,5,6,9,10,15:17)]

for (edit_name in list_to_edit){
  grid.gedit(edit_name, gp = gpar(col="black"))
}


# Manually add axis and label for the dendrogram
decorate_row_dend(" ", {
  vp = current.viewport()
  xscale = vp$xscale
  grid.xaxis(at = xscale[2] - c(0,0.25,0.5,0.75,1)*0.50, label = c(0,0.25,0.5,0.75,1), main=F)
  grid.text("Distance", x = 0.5, y =1.09)
})

decorate_row_dend("  ", {
  vp = current.viewport()
  xscale = vp$xscale
  grid.xaxis(at = xscale[2] - c(0,0.25,0.5,0.75,1)*0.62, label = c(0,0.25,0.5,0.75,1), main=F)
  grid.text("Distance", x = 0.5, y =1.09)
})

decorate_row_dend("   ", {
  vp = current.viewport()
  xscale = vp$xscale
  grid.xaxis(at = xscale[2] - c(0,0.25,0.5,0.75,1)*1.40, label = c(0,0.25,0.5,0.75,1), main=F)
  grid.text("Distance", x = 0.5, y =1.09)
})


# 0.48, 0.64, 1.48



###### Supplementary figure 2 ######

plot_row<- plot_grid(PCA_plot,
                     Table_features,
                     nrow = 1,
                     rel_widths = c(1.5,1))

label_row<- plot_grid(NULL,
                      NULL,
                      nrow = 1,
                      rel_widths = c(1.5,1),
                      labels = c("A","B"),
                      label_size = 20)


plot_grid(label_row,
          plot_row,
          nrow=2,
          rel_heights = c(1,20))



###### Supplementary figure 3 ######

# write the GLUTag rownames in red
cols <- rep("black",dim(heatmap_all@matrix)[1])
cols[which(row.names(heatmap_all@matrix) ==  "Mouse cell line GLUTag")] <- "red"
heatmap_all@row_names_param$gp <- grid::gpar(col=cols, fontsize=8)

cols <- rep("black",dim(heatmap_all@matrix)[2])
cols[which(colnames(heatmap_all@matrix) ==  "Mouse cell line\nGLUTag")] <- "red"
heatmap_all@column_names_param$gp <- grid::gpar(col=cols, fontsize=8)


cols <- rep("black",dim(heatmap_ligandsreceptors@matrix)[1])
cols[which(row.names(heatmap_ligandsreceptors@matrix) ==  "Mouse cell line GLUTag")] <- "red"
heatmap_ligandsreceptors@row_names_param$gp <- grid::gpar(col=cols, fontsize=8)

cols <- rep("black",dim(heatmap_ligandsreceptors@matrix)[2])
cols[which(colnames(heatmap_ligandsreceptors@matrix) ==  "Mouse cell line\nGLUTag")] <- "red"
heatmap_ligandsreceptors@column_names_param$gp <- grid::gpar(col=cols, fontsize=8)


cols <- rep("black",dim(heatmap_gpcrionchannel@matrix)[1])
cols[which(row.names(heatmap_gpcrionchannel@matrix) ==  "Mouse cell line GLUTag")] <- "red"
heatmap_gpcrionchannel@row_names_param$gp <- grid::gpar(col=cols, fontsize=8)

cols <- rep("black",dim(heatmap_gpcrionchannel@matrix)[2])
cols[which(colnames(heatmap_gpcrionchannel@matrix) ==  "Mouse cell line\nGLUTag")] <- "red"
heatmap_gpcrionchannel@column_names_param$gp <- grid::gpar(col=cols, fontsize=8)


# Convert heatmap objects to a format compatible with cowplot
heatmap_all <- grid::grid.grabExpr(draw(heatmap_all))
heatmap_ligandsreceptors <- grid::grid.grabExpr(draw(heatmap_ligandsreceptors))
heatmap_gpcrionchannel <- grid::grid.grabExpr(draw(heatmap_gpcrionchannel))

# Generate plot 
plot_grid(heatmap_all,heatmap_ligandsreceptors,heatmap_gpcrionchannel,nrow=3, labels = c("A","B","C"))



# Only add the lines needed (these are not set beforehand, else black lines would appears for the heatmap even if tits dimension has been set to 0)
list_to_edit <- grid.ls(grid.force(),flatten = T, print=F) %>%
  unlist() %>%
  as.data.frame() %>%
  filter(grepl("GRID.rect",.)) %>%
  pull(.)
list_to_edit <- list_to_edit[-c(1,10,13,22,25,34)]

for (edit_name in list_to_edit){
  grid.gedit(edit_name, gp = gpar(col="black"))
}




###### Supplementary figure 4 ######



# Generate plot 
title_first <- ggdraw() + 
  draw_label(
    "Metabotropic receptor subunits expression in L-cells and related cell models",
    fontface = 'bold',
    hjust = 0.5,
    size = 16
  )

title_second <- ggdraw() + 
  draw_label(
    "AMPA receptor subunits expression in L-cells and related cell models",
    fontface = 'bold',
    hjust = 0.5,
    size = 16
  )

title_third <- ggdraw() + 
  draw_label(
    "Kainate receptor subunits expression in L-cells and related cell models",
    fontface = 'bold',
    hjust = 0.5,
    size = 16
  )

first_row <- plot_grid(list_geneplots[["Metabo_1"]],NULL,list_geneplots[["Metabo_2"]],rel_widths = c(1, 0.05, 1),nrow = 1)
second_row <- plot_grid(list_geneplots[["AMPA_1"]],NULL,list_geneplots[["AMPA_2"]],rel_widths = c(1, 0.05, 1),nrow = 1)
third_row <- plot_grid(list_geneplots[["Kaina_1"]],NULL,list_geneplots[["Kaina_2"]],rel_widths = c(1, 0.05, 1),nrow = 1)


plot_grid(title_first,
          first_row,
          title_second,
          second_row,
          title_third,
          third_row,
          nrow=6,
          rel_heights = c(0.1, 0.9, 0.1, 0.9, 0.1, 0.9),
          labels = c("A","","B","","C",""))











