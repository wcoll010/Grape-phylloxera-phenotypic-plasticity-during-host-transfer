library(readxl)
library(tidyverse)
library(limma)
library(edgeR)
library(statmod)
library(data.table)
library(gplots)
library(dplyr)

setwd("/Users/natecollison/Desktop/DV_RNAseq_10.2023")

getwd()

# read in the gene count matrix generated from STAR alignment
STAR_counts <- read.delim("stargenecount.txt")
count_data_merged <- STAR_counts

# add gene ID as row name
rownames(count_data_merged)<-count_data_merged$gene_id
rownames(count_data_merged)

## remove first row to get only count data
count_data_merged<- count_data_merged[-c(1)]

## check library sizes (raw gene counts)
library.sizes <- colSums(count_data_merged)

library.sizes  
### results:
#G1103_1  G1103_2  G1103_3  G1103_4 G10114_1 G10114_2 G10114_3 G10114_4   G99R_1   G99R_2   G99R_3 
#20830101 20804568 22049301 21922616 23932053 23847363 25969959 18039846 22288110 24885182 22660208 
#G99R_4     G0_1     G0_2     G0_3     G0_4     G1_1     G1_2     G1_3     G1_4    GS4_1    GS4_2 
#21913379 24357207 24423235 15185653 26233968 25882628 26000709 22615726 23773057 22011289 23383086 
#GS4_3    GS4_4 
#18173430 23507863 




range(library.sizes)
# results: 15,185,653 to 26,233,968

# make a bar plot
barplot(library.sizes, las = 2)

# prepare data frame before creating the edgeR object DGEList
targetdata<-data.frame(FileName=colnames(count_data_merged),Group=c(rep("G1103",4), rep("G101_14",4), rep("G99R",4), rep("G0_cab33",4),rep("G1_cab33",4),rep("GS04",4)))
targetdata$Gp2<-factor(targetdata$Group)
targetdata$Gp2<-as.numeric(targetdata$Gp2)
targetdata

x<-DGEList(counts=as.matrix(count_data_merged),lib.size=library.sizes,group=targetdata$Group)


class(x) 
#[1] "DGEList"
#attr(,"package")
#[1] "edgeR"

dim(x) 
#25814    24

colnames(x)

samplenames <- substring(colnames(x), 1, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
colnames(x)
rownames(x)

group <- as.factor(c( "G1103", "G1103", "G1103", "G1103","G101_14", "G101_14", "G101_14", "G101_14",  "G99R", "G99R", "G99R", "G99R","G0_cab33", "G0_cab33", "G0_cab33", "G0_cab33","G1_cab33", "G1_cab33", "G1_cab33", "G1_cab33", "GS04", "GS04","GS04", "GS04"))


x$samples$group <- group
x$samples

## filter out gene with insufficient counts using the filterByExpr function
?filterByExpr

keep.exprs <- filterByExpr(x, group=group) 
x <- x[keep.exprs,, keep.lib.sizes=FALSE] 
dim(x) 
# 13,125    24


# make multi-dimentional scaling plot (MDS) using TMM-normalized gene count data
?calcNormFactors
x <- calcNormFactors(x, method = "TMM")
plotMDS(x)

x$samples$norm.factors

# load libraries
library(RColorBrewer)
library(vegan) # for adding circles to the mds plot

# Make color palette without yellow
my_colors <- brewer.pal(8,"Set2")
my_colors <- my_colors[-5]
my_colors <- c(my_colors[length(my_colors)], my_colors[-length(my_colors)])

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,1))
col.group <- group

levels(col.group) <- my_colors
col.group <- as.character(col.group)

# Create MDS plot with points only and store the result
mds <- plotMDS(lcpm, col=col.group, pch=20) # pch=16 gives solid circles

# Extract MDS coordinates (x and y dimensions)
mds_coords <- cbind(mds$x, mds$y)

# Add ellipses around each treatment group
ordiellipse(mds_coords, group, col=my_colors, kind="ehull", draw="lines", lwd=2)

# Add a legend
legend("topright", legend=levels(group), col=my_colors, pch=16, cex=0.8)



### Differential expression (DEG) analysis
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

# make contrasts
#compare all hosts
contr.matrix <- makeContrasts(
  G1_cab33vsG101_14 = G1_cab33 - G101_14,
  G1_cab33vsG99R = G1_cab33 - G99R,
  G1_cab33vsG1103 = G1_cab33 - G1103,
  G1_cab33vsGS04 = G1_cab33 - GS04,
  G1_cab33vsG0_cab33 = G1_cab33 - G0_cab33,
  G101_14vsG99R = G101_14 - G99R,
  G101_14vsG1103 = G101_14 - G1103,
  G101_14vsGS04 = G101_14 - GS04,
  G99RvsG1103 = G99R - G1103,
  G99RvsGS04 = G99R - GS04,
  G1103vsGS04 = G1103 - GS04,
  G0_cab33vsG101_14 = G0_cab33 - G101_14,
  G0_cab33vsG99R = G0_cab33 - G99R,
  G0_cab33vsG1103 = G0_cab33 - G1103,
  G0_cab33vsGS04 = G0_cab33 - GS04,
  levels = colnames(design)
)

contr.matrix

# run voom model fitting
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design) 
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) 
efit <- eBayes(vfit)
plotSA(efit)
# table of DEG
summary(decideTests(efit)) 


# add fold-change filter
?treat
tfit <- treat(vfit, fc=1.5) 
dt <- decideTests(tfit) 
summary(dt)


# plot mean-difference MD
?plotMD
par(mfrow=c(1,1))
plotMD(tfit, column=1, status=dt[,1], main="", xlim=c(-8,13))

## plot a venn diagram

vennDiagram(dt[, c(4, 8, 10, 11)], 
            circle.col = c("turquoise", "salmon", "blue4", "green3"), 
            category.names = colnames(dt)[c(4, 8, 10, 11)],
            cex = 1,
            cex.category = 0.5) 



## annotate DEG list by merging with gene annotation file (gene_annotations.csv, transdata2.csv)
tpmdata<-read.csv("gene_annotations.csv", header=TRUE,stringsAsFactors=T)
tpmdata<-tpmdata[,c(1:15)]
dt_df <- as.data.frame(dt)

dt_df



### Find gene that are differentially expression across different host genotypes, i.e. host-responsive genes

# Step 1: Keep only genes that are differentially expressed (-1 or 1 in at least one column)
diff_expressed_df <- dt_df %>%
  filter(if_any(everything(), ~ . == -1 | . == 1))

#remove genes that are DE in GOvsG1 contrast
filtered_DEG_df <- diff_expressed_df %>%
  filter(G1_cab33vsG0_cab33 == 0)

# add geneID (fasta3.1ID) as a column in the dataframe
filtered_DEG_df$fasta3.1ID <- rownames(filtered_DEG_df)

# add gene annotation info via merging with tpmdata object (annotation file)
data_frame_merge <- merge(tpmdata,filtered_DEG_df, by = 'fasta3.1ID', all = FALSE)

#write data frame as excel file
writexl::write_xlsx(data_frame_merge, "host_responsive_DEGs_minusG0_03182025.xlsx")

# write a dataframe for the genes that are DE in G0 vs G1 
DEG_G1_vs_G0 <- diff_expressed_df %>%
  filter(G1_cab33vsG0_cab33 != 0)
DEG_G1_vs_G0$fasta3.1ID <- rownames(DEG_G1_vs_G0)
data_frame_merge_2 <- merge(tpmdata,DEG_G1_vs_G0, by = 'fasta3.1ID', all = FALSE)
writexl::write_xlsx(data_frame_merge_2, "DEGs_G1_minus_G0_03192025.xlsx")




# write a new data frame removing genes that are DE in G0vsG1 , but all other genes are retained in the data
filtered_DEG_df_2 <- dt_df %>%
  filter(G1_cab33vsG0_cab33 == 0)
# 11,753

## plot a venn diagram
# for G1_cab33 contrasts
vennDiagram(filtered_DEG_df_2[, 1:4], 
            circle.col = c("turquoise", "salmon", "blue4", "green3"), 
            category.names = colnames(filtered_DEG_df_2)[1:4],
            cex = 1,
            cex.category = 0.5) 

# 3,432 genes are DEG out of 3706 total, 92.61% 

# for S04 host genotype
vennDiagram(filtered_DEG_df_2[, c(4, 8, 10, 11)], 
            circle.col = c("turquoise", "salmon", "blue4", "green3"), 
            category.names = colnames(filtered_DEG_df_2)[c(4, 8, 10, 11)],
            cex = 1,
            cex.category = 0.5) 


# write a data frame for genes that are upregulated on Cab33 G1
df <- read_xlsx("host_responsive_DEGs_minusG0_03182025.xlsx")

# Filter for rows where any of columns 16:19 have a value of 1 (upregulated)
df_upregulated <- df[apply(df[, 16:19] == 1, 1, any), ]
#2,667 genes

# Save to file
writexl::write_xlsx(df_upregulated, "genes_upregulated_on_G1cab33.xlsx")


# write data frame for genes that are up-regulated or down regulated on Cab33 G1
df_de_cab33 <- df[apply(df[, 16:19], 1, function(x) any(x %in% c(1, -1))), ]
# 3,432 genes

writexl::write_xlsx(df_de_cab33, "genes_DE_on_G1cab33.xlsx")



### making pie charts for gene functional categories

library(tidyverse)
library(readxl)

# read in files
data_frame_merge <-read_xlsx("host_responsive_DEGs_minusG0_03182025.xlsx")
tpmdata<-read.csv("transdata2.csv", header=TRUE,stringsAsFactors=T)
tpmdata<-tpmdata[,c(1:15)]

samplesize <- tpmdata %>% 
  group_by(AnnotGroup)%>% 
  tally()
samplesize

samplesize <- data_frame_merge_2 %>% 
  group_by(AnnotGroup)%>% 
  tally()
samplesize

#2 "Circadian rhythm"                                                 16
#3 "Cuticle and chitin related proteins"                             116
#4 "Detoxification of plant and other environmental molecules"       222
#5 "Developmental genes"                                             135
#6 "Effectors"                                                      2743
#7 "Epigenetic machinery"                                              5
#8 "General metabolism"                                               14
#9 "Immune system"                                                     2
#10 "null"                                                          20998
#11 "Odorant or gustatory receptors and others associated proteins"   107
#12 "Other"                                                          1280
#13 "Ribosomal proteins"                                              155
#14 "Small RNA machinery"                                              12
#15 "Transporters"                                                      3


## make pie chart
library(graphics)
library(RColorBrewer)

## switch this for plotting other datasets
data<-data_frame_merge
data<-tpmdata

data$AnnotGroup <- as.character(data$AnnotGroup)
data$AnnotGroup <- sub("Detoxification of plant and other environmental molecules", "Detoxification", data$AnnotGroup)
data$AnnotGroup <- sub("Odorant or gustatory receptors and others associated proteins", "Odorant or gustatory perception", data$AnnotGroup)

counts <- table(data$AnnotGroup)

# Calculate percentages
percentages <- prop.table(counts) * 100

# Define the desired order of levels
desired_order <- c("null","Detoxification", "Other", "Effectors","Developmental genes","Cuticle and chitin related proteins", "Odorant or gustatory perception","Ribosomal proteins", "Circadian rhythm","Transporters","Small RNA machinery","Epigenetic machinery","General metabolism","Immune system")

# Reorder the factor levels based on the desired order
data$AnnotGroup <- factor(data$AnnotGroup, levels = desired_order)

# Choose a different color palette from RColorBrewer
colors <- brewer.pal(length(desired_order), "Set3")

# Create a table of counts
counts <- table(data$AnnotGroup)

# Create a pie chart with the selected color palette
par(mfrow=c(1,1))
par(cex = 1, mar = c(1, 1, 1, 1)) # Adjust the margin size

pie(counts, labels = NA, col = colors, main = NA, cex.main = 0.8)

pie(counts, labels = ifelse(counts / sum(counts) * 100 > 0.2, names(counts), ""), col = colors, main = "Distribution of major annotation groups: host responsive genes (fc 1.5) ", cex.main = 0.8)
legend("topleft", legend = c("null", "Detoxification", "Other", "Predicted effectors", "Developmental genes","Cuticle and chitin"), fill = colors, title = "Annotation Groups")



## Perform Fisher's exact test for annotation group enrichment 

df_genes <- data_frame_merge
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

fischer_test <- as.data.frame(fisher_test_result)

# using BH adjusted p values 

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

# Print the results
print(results)
results <- as.data.frame(results)
writexl::write_xlsx(results, "fisher_test_BHad_1.5fc_DEGs_minusG0.xlsx")

##
