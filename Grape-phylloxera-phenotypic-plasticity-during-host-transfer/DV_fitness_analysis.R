setwd("/Users/natecollison/Desktop/Nabity Lab")

getwd()


# Load necessary library
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(dplyr)

# load data
data <- read_xlsx("June performance data.xlsx")

# Replace genotype names using mutate and recode
data <- data %>%
  mutate(genotype = recode(genotype,
                           "Cab 33" = "G1_cab33",
                           "101-14" = "G101_14",
                           "1103P" = "G1103",
                           "99R" = "G99R",
                           "SO4" = "GS04"))

# Verify changes
head(data)

# Make color palette from Set1 without yellow
my_colors <- brewer.pal(8,"Set2")
my_colors <- my_colors[-5]  # Remove the yellow color

# Create the violin plot with updated genotype names
ggplot(data, aes(x = genotype, y = eggs, fill = genotype)) +
  geom_violin(trim = FALSE) +
  #geom_jitter(width = 0.2, size = 1.5, color = "black") +
  labs(x = "host genotype", y = "eggs per gall") +
  theme_minimal() +
  scale_fill_manual(values = my_colors) +
  theme(legend.title = element_blank())

?theme
##
# Create the box plot with custom colors
ggplot(data, aes(x = genotype, y = eggs, fill = genotype)) +
  geom_boxplot() +  # Box plot
  geom_point(position = position_jitter(width = 0.2), shape = 1, size = 2) +  # Add hollow data points
  labs(x = "Genotype", y = "Egg count", fill = NULL) +  # Remove legend title
  theme_minimal() +  # Minimal theme for a clean plot
  scale_fill_manual(values = my_colors) +  # Apply your custom color palette
  theme(legend.title = element_blank())  # Remove the legend title



## statistical tests
# Check normality for each genotype
shapiro_test <- data %>% 
  group_by(genotype) %>% 
  summarise(p_value = shapiro.test(eggs)$p.value)

print(shapiro_test)

# If data is normally distributed for each genotype (p > 0.05), run ANOVA
# data not normally dist for Cab33 and GS04
# A tibble: 5 × 2
#genotype    p_value
#<chr>         <dbl>
#  1 G101_14  0.910     
#2 G1103    0.0801    
#3 G1_cab33 0.00000598
#4 G99R     0.128     
#5 GS04     0.0000450 

# If the data is not normally distributed (p < 0.05), run Kruskal-Wallis test
kruskal_result <- kruskal.test(eggs ~ genotype, data = data)
print(kruskal_result)
#Kruskal-Wallis rank sum test

#data:  eggs by genotype
#Kruskal-Wallis chi-squared = 175.51, df = 4, p-value <  2.2e-16


# If Kruskal-Wallis is significant, perform pairwise Wilcoxon tests
pairwise_wilcox <- pairwise.wilcox.test(data$eggs, data$genotype, p.adjust.method = "BH")
print(pairwise_wilcox)

# reprint box plot with stats

library(ggsignif)

ggplot(data, aes(x = genotype, y = eggs, fill = genotype)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), shape = 1, size = 2) +
  labs(x = "Genotype", y = "Egg count", fill = NULL) +
  theme_minimal() +
  scale_fill_manual(values = my_colors) +
  theme(legend.title = element_blank()) +
  
  # Add significance annotations between groups
  geom_signif(comparisons = list(c("G1_cab33", "G101_14"), 
                                 c("G1_cab33", "G1103"),
                                 c("G1_cab33", "G99R"),
                                 c("G1_cab33", "GS04"),
                                 c("G101_14", "GS04"),
                                 c("G1103", "GS04"),
                                 c("G99R", "GS04")),
              map_signif_level = TRUE)



##