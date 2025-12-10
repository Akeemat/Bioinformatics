getwd()
list.files()

#Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

#read columns as character to prevent parsing errors
gwas_covid_data <- read_tsv("gwas-association-downloaded_2025-11-15-MONDO_0100096.tsv", col_types = cols(.default = "c"))

#check data
head(gwas_covid_data)
str(gwas_covid_data)
colnames(gwas_covid_data)
nrow(gwas_covid_data)
#convert numeric columns to numeric
numeric_columns <- c(
  "P-VALUE",
  "PVALUE_MLOG",
  "CHR_POS",
  "MERGED",
  "SNP_ID_CURRENT",
  "UPSTREAM_GENE_DISTANCE",
  "DOWNSTREAM_GENE_DISTANCE",
  "INTERGENIC",
  "OR or BETA"
)

gwas_covid_data[numeric_columns] <- lapply(gwas_covid_data[numeric_columns], \(x) suppressWarnings(as.numeric(x)))

# Filter for all relevant COVID-19 traits
relevant_traits <- c(
  # General COVID-19
  "COVID-19",
  "SARS-CoV-2 infection",
  "COVID-19 susceptibility",
  
  # Susceptibility
  "COVID-19 (covid vs negative)",
  "Resistance to COVID-19 infection (Exposed negative vs positive)",
  
  # Hospitalization
  "COVID-19 (hospitalized covid vs population)",
  "COVID-19 (hospitalized vs not hospitalized)",
  "COVID-19 (hospitalized vs population)",
  "COVID-19 (hospitalized vs tested, not hospitalized)",
  
  # Severity - General
  "COVID-19 (severe vs non-severe with positive test)",
  "COVID-19 (severe vs population)",
  "COVID-19 (severe vs tested, not severe)",
  "Severe COVID-19 infection",
  "Severe COVID-19 disease",
  
  # Severity - Specific outcomes
  "COVID-19 (critical illness vs population)",
  "COVID-19 (covid pneumonia vs population)",
  "COVID-19 (covid respiratory support vs population)",
  "COVID-19 (severe respiratory symptoms vs population)",
  "COVID-19 with respiratory failure",
  "Severe COVID-19 infection with respiratory failure (analysis I)",
  
  # Mortality
  "COVID-19 death (death vs population)",
  "COVID-19 death (death vs tested and survived)",
  "Mortality in COVID-19",
  
  # Symptom presentation
  "Asymptomatic COVID-19 infection",
  "COVID-19 (recovered vs asymptomatic)",
  "COVID-19 (symptomatic vs paucisymptomatic)",
  
  # Other relevant
  "Morbidity in COVID-19 (deceased vs asymptomatic)"
)

gwas_filtered <- gwas_covid_data %>%
  filter(`DISEASE/TRAIT` %in% relevant_traits)

#Filter further for significant associations
gwas_significant <- gwas_filtered %>%
  filter(`P-VALUE` < 5e-8)
nrow(gwas_significant)

#Generate summary of significant associations by chromosome
chromosome_summary <- gwas_significant %>%
  filter(!is.na(CHR_ID), CHR_ID != "") %>%
  group_by(CHR_ID) %>%
  summarise(
    SNP_count = n(),
    Min_P_value = min(`P-VALUE`, na.rm = TRUE),
    Max_OR_BETA = max(`OR or BETA`, na.rm = TRUE),
    Top_genes = paste(unique(MAPPED_GENE[!is.na(MAPPED_GENE)]), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(SNP_count))

#Create a summary table of significant associations by gene

gene_summary <- gwas_significant %>%
  filter(!is.na(MAPPED_GENE), MAPPED_GENE != "") %>%
  mutate(MAPPED_GENE = str_replace_all(MAPPED_GENE, " - ", ",")) %>%
  separate_rows(MAPPED_GENE, sep = ",") %>%
  group_by(MAPPED_GENE) %>%
  summarise(
    SNP_Count = n(),
    Chromosomes = paste(unique(CHR_ID), collapse = ", "),
    Min_P_value = if(all(is.na(`P-VALUE`))) NA else min(`P-VALUE`, na.rm = TRUE),
    Max_OR_BETA = if(all(is.na(`OR or BETA`))) NA else max(`OR or BETA`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(SNP_Count))

head(gene_summary)

#Generate barplot of significant associations by gene

library(ggplot2)

# Select top 12 genes
top12_genes <- gene_summary %>%
  slice_max(SNP_Count, n = 12)

ggplot(top12_genes,
       aes(x = reorder(MAPPED_GENE, SNP_Count), y = SNP_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = SNP_Count),
            hjust = -0.1, size = 4) +
  coord_flip() +
  labs(
    x = "Gene",
    y = "Number of Significant SNPs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    axis.text.x  = element_text(color = "black"),
    axis.text.y  = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black")
  ) +
  ylim(0, max(top12_genes$SNP_Count) * 1.1)

# Further analysis for signifant SNPs on ACE2
# Filter rows where ACE2 is the mapped gene
ace2_snps_significant <- gwas_significant %>%
  filter(str_detect(MAPPED_GENE, "ACE2"))

#Export cleaned data into csv 
getwd()
write.csv(ace2_snps_significant, "ace2_snps_significant.csv", row.names = FALSE)
write.csv(gwas_significant, "gwas_significant.csv", row.names = FALSE)
