#Set and confirm working directory
getwd()
setwd("C:/Users/User/Documents/BMI_5330") 
list.files()

library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(tidyverse)

# Read in file
dbvar_entries <- read_file("dbvar_result.txt")

#confirm file details
head(dbvar_entries)
str(dbvar_entries)
nrow(dbvar_entries) 

# Split into blocks using the numbering pattern
blocks <- unlist(str_split(dbvar_entries, "\\n(?=\\d+\\. )"))

# Remove empty blocks and blocks with no numbers
blocks <- blocks[nchar(blocks) > 0 & str_detect(blocks, "^\\d+\\.")]

# Create function to extract each entry
extract_field <- function(text, pattern) {
  m <- str_match(text, pattern)
  if (is.na(m[1,2])) return(NA)
  return(str_trim(m[1,2]))
}

# Parse each block
extracted_variants <- lapply(blocks, function(b) {
  list(
    entry_number     = extract_field(b, "^(\\d+)\\."),
    variant_id       = extract_field(b, "(nsv\\d+|esv\\d+)"),
    variant_type     = extract_field(b, "Variant type:\\s*(.*)"),
    associated_study = extract_field(b, "Associated study:\\s*(.*)"),
    organism         = extract_field(b, "Organism:\\s*(.*)"),
    genes            = extract_field(b, "Gene\\(s\\) in region:\\s*(.*)"),
    location_id      = extract_field(b, "ID:\\s*(\\d+)")
  )
})

variant_table <- bind_rows(extracted_variants)

# Convert numeric columns
variant_table <- variant_table %>%
  mutate(entry_number = as.integer(entry_number),
         location_id  = as.integer(location_id))

#Deduplicate and filter for only variations related to ACE2 and human

dbvar_table_unique <- variant_table %>%
  distinct() %>%                                    # remove exact duplicates
  filter(str_detect(genes, "ACE2")) %>%            # keep only rows with 'ACE2'
  filter(str_to_lower(organism) == "human")       # keep only human entries

# Quick check
nrow(dbvar_table_unique)
head(dbvar_table_unique)

#Determine frequency of variant types
table(dbvar_table_unique$variant_type)
prop.table(table(dbvar_table_unique$variant_type))

#Summarize mutation types
dbvar_varianttype_summary <- as.data.frame(table(dbvar_table_unique$variant_type))
colnames(dbvar_varianttype_summary) <- c("Variant_Type", "Frequency")
dbvar_varianttype_summary <- dbvar_varianttype_summary[order(dbvar_varianttype_summary$Frequency, decreasing = TRUE), ]
sum(dbvar_varianttype_summary$Frequency)

#Make a barplot of structural variations
library (ggplot2)
library(forcats)

# Reorder variant types ascending 
dbvar_varianttype_summary$Variant_Type <- fct_reorder(
  dbvar_varianttype_summary$Variant_Type,
  dbvar_varianttype_summary$Frequency,
  .desc = FALSE
)

ggplot(dbvar_varianttype_summary,
       aes(x = Variant_Type, y = Frequency)) +
  geom_col(fill = "#377eb8") +  # nice blue
  geom_text(aes(label = Frequency),
            hjust = -0.2, size = 4) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Variant type",
    y = "Number of variants"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(color = "black"),
    axis.text.y     = element_text(color = "black"),
    axis.title.x    = element_text(color = "black"),
    axis.title.y    = element_text(color = "black"),
    plot.title      = element_text(color = "black", face = "bold")
  )

#Export cleaned data into csv 
write.csv(dbvar_table_unique, "dbvar_table_unique.csv", row.names = FALSE)
