# load packages/libraries
library(clusterProfiler)
library(enrichplot)
library(data.table)
library(tidyr)
library(readxl)
library(tidyverse)

# set working directory
setwd("/Users/natecollison/Desktop/DV_RNAseq_10.2023")
getwd()



### Perform KEGG term enrichment analysis and make dot plots - what KEGG terms are enriched among host-responsive genes

# load KEGG annotation data
eggnog_data <- read_xlsx("MM_9qtlz_3b.emapper.annotations.xlsx")

# remove extra flags on gene ID (-PA) so gene IDs match the gene expression data
eggnog_data$query <- gsub("-PA","",as.character(eggnog_data$query))
colnames(eggnog_data)[1] ="ID"


# load the list of host-responsive genes
DEG_list <- read_xlsx("host_responsive_genes_fc1.5_09132024.xlsx") 
DEG_list<- DEG_list[-c(1)]

# make a list genes down-regulated
DOWN_DEGs <- subset(DEG_list, DEG_list[, 15] == -1 | DEG_list[, 16] == -1 | DEG_list[, 17] == -1 | DEG_list[, 18] == -1)
# down degs = 884 genes

# up-regulated genes
UP_DEGs <- subset(DEG_list, DEG_list[, 15] == 1 | DEG_list[, 16] == 1 | DEG_list[, 17] == 1 | DEG_list[, 18] == 1)
# UP degs = 3,933 genes


# switch DEG_list_KO_merged for UP and DOWN DEGs sets
DEG_list_KO_merged <- merge(DOWN_DEGs,eggnog_data, by = 'ID', all = FALSE)
DEG_KOs <- DEG_list_KO_merged[c(1,35)] # get only gene name and KEGG term for new dataframe

# get gene vector for DEG KOs
DEG_KO_vector <- DEG_KOs$ID

# create KEGG background object (all genes and KEGG terms, ie KEGG universe)
KO_universe <- eggnog_data[c(1,12)]

KO_universe$KEGG_ko <- gsub( "ko:", "", as.character(KO_universe$KEGG_ko))

# expand data frames so each KO term in a unique row
KO_universe <- separate_rows(KO_universe, KEGG_ko, sep = ",")
KO_universe <- KO_universe[,c(2,1)]


# run clusterProfiler KEGG enrichment
KO_enrichment_results <- enricher(DEG_KO_vector, TERM2GENE=KO_universe, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)

?enricher

# make a simple dot plot of the results and save it as a .tiff file
pdf(file="KO_enrich_DOWN_DEGvsCab33_1.5fc_dotplot.pdf",width=6,height=4)
dotplot(KO_enrichment_results, showCategory=100)
dev.off()






### KEGG term enrichment analysis for gene co expression modules 

setwd("/Users/natecollison/Desktop/DV_RNAseq_10.2023")
getwd()


## load KEGG term list / KO universe or genome background
eggnog_data <- read_xlsx("MM_9qtlz_3b.emapper.annotations.xlsx")
eggnog_data$query <- gsub("-PA","",as.character(eggnog_data$query))
colnames(eggnog_data)[1] ="ID"

gene_mod_assignments <- read_xlsx("gene_mod_assign_signed_16thpwr_unmerged_07162024.xlsx")

darkorange2_df <- gene_mod_assignments[gene_mod_assignments$ModuleAssignment == "darkorange2", ] # change for diff modules
dim(darkorange2_df)

thistle_df <- gene_mod_assignments[gene_mod_assignments$ModuleAssignment == "thistle", ]
dim(thistle_df)

darkolivegreen4_df <- gene_mod_assignments[gene_mod_assignments$ModuleAssignment == "darkolivegreen4", ]
dim(darkolivegreen4_df)

antiquewhite4_df <- gene_mod_assignments[gene_mod_assignments$ModuleAssignment == "antiquewhite4", ]
dim(antiquewhite4_df )

##
gene_list_KO_merged <- merge(antiquewhite4_df,eggnog_data, by = 'ID', all = FALSE) # change the dataframe here (e.g. antiquewhite4_df) to perfrom analysis on other modules)
gene_KOs <- gene_list_KO_merged[c(1,27)] # get only gene name and KEGG term for new dataframe

KO_universe <- eggnog_data[c(1,12)]

KO_universe$KEGG_ko <- gsub( "ko:", "", as.character(KO_universe$KEGG_ko))

# expand data frames so each KO term in a unique row
KO_universe <- separate_rows(KO_universe, KEGG_ko, sep = ",")
KO_universe <- KO_universe[,c(2,1)]

# get gene vector for DEG KOs
gene_KO_vector <- gene_KOs$ID

# run clusterProfiler KEGG term enrichment
?enricher
KO_enrichment_results <- enricher(gene_KO_vector, TERM2GENE=KO_universe, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)
dotplot(KO_enrichment_results, showCategory=100)

# make a simple dot plot of the results and save it as a pdf file
pdf(file="KO_enrich_MEantiquewhite4.pdf",width=6,height=4)
dotplot(KO_enrichment_results, showCategory=100)
dev.off()






## Perform Fisher's exact test for functional group enrichment of genes up-regulated and down-regulated on non-native host (Cab33 minus other hosts(pairwise))

#load gene annotation file 
tpmdata<-read.csv("transdata2.csv", header=TRUE,stringsAsFactors=T)
tpmdata<-tpmdata[,c(1:15)]

DEG_list <- read_xlsx("common_DEGS_Cab33vsall.xlsx") 
DEG_list<- DEG_list[-c(1)]


# make a list genes down-regulated
DOWN_DEGs <- subset(DEG_list, DEG_list[, 15] == -1 | DEG_list[, 16] == -1 | DEG_list[, 17] == -1 | DEG_list[, 18] == -1)

# up-regulated genes
UP_DEGs <- subset(DEG_list, DEG_list[, 15] == 1 | DEG_list[, 16] == 1 | DEG_list[, 17] == 1 | DEG_list[, 18] == 1)

# check number of genes per functional group
samplesize <- tpmdata %>% 
  group_by(AnnotGroup)%>% 
  tally()
samplesize

samplesize <- UP_DEGs %>% 
  group_by(AnnotGroup)%>% 
  tally()
samplesize

samplesize <- DOWN_DEGs %>% 
  group_by(AnnotGroup)%>% 
  tally()
samplesize

# write R objects to new dataframe
df_genes <- UP_DEGs
df_genes <- DOWN_DEGs
df_genes <- DEG_list

df_all_genes <- tpmdata

# Get unique annotation groups
unique_groups <- unique(df_all_genes$AnnotGroup)

# Loop through each module and annotation group
for (group in unique_groups) {
  # Count genes in the module belonging to the current annotation group
  group_count_module <- sum(df_genes$AnnotGroup == group)
  
  # Count genes in the entire dataset belonging to the current annotation group
  group_count_all <- sum(df_all_genes$AnnotGroup == group)
  
  # Create a contingency table
  contingency_table <- matrix(c(group_count_module, sum(df_genes$AnnotGroup != group),
                                group_count_all - group_count_module, sum(df_all_genes$AnnotGroup != group)),
                              nrow = 2, byrow = TRUE)
  
  # Perform Fisher's exact test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Print results
  cat("Annotation Group:", group, "\n")
  print(fisher_test_result)
  cat("\n")
}


#####
# using BH-adjusted p values 
df_genes <- UP_DEGs
df_genes <- DOWN_DEGs
df_genes <- DEG_list

df_all_genes <- tpmdata

# Get unique annotation groups
unique_groups <- unique(df_all_genes$AnnotGroup)

# Initialize vectors to store p-values and odds ratios
p_values <- numeric(length(unique_groups))
odds_ratios <- numeric(length(unique_groups))

# Perform Fisher's exact test for each group and store the p-values and odds ratios
for (i in seq_along(unique_groups)) {
  group <- unique_groups[i]
  
  # Count genes in the module belonging to the current annotation group
  group_count_module <- sum(df_genes$AnnotGroup == group)
  
  # Count genes in the entire dataset belonging to the current annotation group
  group_count_all <- sum(df_all_genes$AnnotGroup == group)
  
  # Create a contingency table
  contingency_table <- matrix(c(group_count_module, sum(df_genes$AnnotGroup != group),
                                group_count_all - group_count_module, sum(df_all_genes$AnnotGroup != group)),
                              nrow = 2, byrow = TRUE)
  
  # Perform Fisher's exact test
  ?fisher.test
  fisher_test_result <- fisher.test(contingency_table)
  
  # Store the p-value and odds ratio
  p_values[i] <- fisher_test_result$p.value
  odds_ratios[i] <- fisher_test_result$estimate
  
  # Print results
  cat("Annotation Group:", group, "\n")
  cat("Raw p-value:", format(fisher_test_result$p.value, scientific = TRUE), "\n")
  cat("95 percent confidence interval:", fisher_test_result$conf.int, "\n")
  cat("odds ratio:", fisher_test_result$estimate, "\n\n")
}

# Adjust the p-values using the Benjamini-Hochberg method
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Create a data frame to store annotation groups, adjusted p-values, and odds ratios
results <- data.frame(
  AnnotationGroup = unique_groups,
  AdjustedPValue = adjusted_p_values,
  OddsRatio = odds_ratios,
  stringsAsFactors = FALSE
)

# View the results
print(results)

# Write results to excel file 
results <- as.data.frame(results)
writexl::write_xlsx(results, "fisher_test_BHad_1.5fc_DEGs.xlsx")



##
