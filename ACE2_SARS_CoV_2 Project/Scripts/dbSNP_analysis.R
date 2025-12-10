#Set and confirm working directory
getwd()
setwd("C:/Users/User/Documents/BMI_5330") 
list.files()
#read into R
dbSNP_variants <- read.delim("snp_result.txt")

#confirm file details
head(dbSNP_variants)
str(dbSNP_variants)
nrow(dbSNP_variants) 

#filter to exclude irrelevant entries
ace2_variants <- dbSNP_variants[dbSNP_variants$X.chr == "X" & grepl("ACE2", dbSNP_variants$gene),]
nrow(ace2_variants)

##check for duplicates
#Count unique IDs
length(unique(ace2_variants$snp_id))

#deduplicate using frequency column
library(tidyverse)
ace2_variants_unique <- ace2_variants %>%
  arrange(snp_id, desc(frequency != "")) %>%  # prioritize non-empty frequency
  distinct(snp_id, .keep_all = TRUE)

#determine frequency of genes 
table(ace2_variants_unique$gene)
prop.table(table(ace2_variants_unique$gene))

#determine frequency of variant types
table(ace2_variants_unique$variant_type)
prop.table(table(ace2_variants_unique$variant_type))

#Summarize mutation types
varianttype_summary <- as.data.frame(table(ace2_variants_unique$variant_type))
colnames(varianttype_summary) <- c("Variant_Type", "Frequency")
sum(varianttype_summary$Frequency)

#Make a barplot of SNPs
library (ggplot2)
varianttype_summary$Variant_Type <- factor(
  varianttype_summary$Variant_Type,
  levels = c("snv", "delins", "del", "ins", "mnv")
)

ggplot(varianttype_summary,
       aes(x = Variant_Type, y = Frequency, fill = Variant_Type)) +
  geom_col() +
  geom_text(aes(label = Frequency),
            vjust = -0.3, size = 4) +
  labs(
    x = "Variant Type",
    y = "Number of variants"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "none",
    axis.text.x      = element_text(color = "black"),
    axis.text.y      = element_text(color = "black"),
    axis.title.x     = element_text(color = "black"),
    axis.title.y     = element_text(color = "black"),
    plot.title       = element_text(color = "black", face = "bold")
  )
       

#Export cleaned data into csv 
getwd()
write.csv(ace2_variants_unique, "ace2_variants_unique.csv", row.names = FALSE)




