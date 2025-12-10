#### GWAS significant variants ####
#Read in extracted variant file
ace2_4_snps <- read.delim("ace2_4snps_by_position.txt", header = FALSE, stringsAsFactors = FALSE)
#Add column names
colnames(ace2_4_snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
View(ace2_4_snps)
#Add snpID column
ace2_4_snps$RSID <- c(
  "rs4830964",  # row 1 (POS = 15515319)
  "rs35697037",  # row 2 (POS = 15542116)
  "rs60768809",  # row 3 (POS = 15547161)
  "rs190509934"   # row 4 (POS = 15620340)
)

#Split list and extract relevant columns
info_list <- strsplit(ace2_4_snps$INFO, ";") #Split the INFO field
extract_info_field <- function(info_list, key) {
  sapply(info_list, function(x) {
    # find the piece that starts with key=
    m <- grep(paste0("^", key, "="), x, value = TRUE)
    if (length(m) == 0) return(NA)
    # remove "key=" and keep just the value
    sub("^[^=]+=", "", m[1])
  }) 
}
fields <- c("AF", "AFR_AF", "EUR_AF", "EAS_AF", "SAS_AF", "AMR_AF")

for (f in fields) {
  ace2_4_snps[[f]] <- as.numeric(extract_info_field(info_list, f))
}
ace2_4_table <- ace2_4_snps[, c("RSID", "CHROM", "POS", "REF", "ALT", "AF", "AFR_AF", "EUR_AF", "EAS_AF", "SAS_AF", "AMR_AF")]
View(ace2_4_table)

#Export cleaned data into csv 
write.csv(ace2_4_table, "ace2_4_table.csv", row.names = FALSE)


#### Variants within ACE2 region ####
#Read in extracted variant file
ace2_region <- read.delim("ace2_region_variants.txt", header = FALSE, stringsAsFactors = FALSE)

colnames(ace2_region) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

# Quick QC check
nrow(ace2_region)
head(ace2_region$INFO)

# Split INFO at ";" into a list
info_list_region <- strsplit(ace2_region$INFO, ";", fixed = TRUE)

#Create function to extract specific fields
extract_info_field <- function(info_list, key) {
  sapply(info_list, function(x) {
    # find the piece that starts with key=
    m <- grep(paste0("^", key, "="), x, value = TRUE)
    if (length(m) == 0) return(NA)        
    sub("^[^=]+=", "", m[1])              
  })
}

#Extract AF and population-specific AFs
fields <- c("AF", "AFR_AF", "EUR_AF", "EAS_AF", "SAS_AF", "AMR_AF")

for (f in fields) {
  ace2_region[[f]] <- as.numeric(extract_info_field(info_list_region, f))
}

#QC check
head(ace2_region[, c("CHROM", "POS", "REF", "ALT",
                     "AF", "AFR_AF", "EUR_AF", "EAS_AF", "SAS_AF", "AMR_AF")])

summary(ace2_region$AF)
summary(ace2_region$AFR_AF)

#Pivot to long format 
AF_pop_long <- ace2_region %>%
  select(AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF) %>%   # just reordered for clarity
  pivot_longer(
    cols = everything(),
    names_to = "Population",
    values_to = "AF"
  ) %>%
  filter(!is.na(AF)) %>%
  
  #change labels
  mutate(
    Population = factor(
      Population,
      levels = c("AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF"),
      labels = c("AFR", "AMR", "EAS", "EUR", "SAS")
    )
  )
#Calculate mean allele frequencies 
AF_means <- AF_pop_long %>%
  group_by(Population) %>%
  summarise(mean_AF = mean(AF, na.rm = TRUE))


#Create a table with mean, min, max values, no of variants per category

library(dplyr)

AF_stats <- AF_pop_long %>%
  group_by(Population) %>%
  summarise(
    n_variants = sum(!is.na(AF)),
    mean_AF    = mean(AF, na.rm = TRUE),
    min_AF     = min(AF, na.rm = TRUE),
    max_AF     = max(AF, na.rm = TRUE)
  ) %>%
  arrange(Population)

#Create barplot comparing population-specific AFs

library(dplyr)
library(ggplot2)

# Assuming AF_pop_long has columns: Population (AFR/AMR/EAS/EUR/SAS) and AF (numeric)
AF_summary <- AF_pop_long %>%
  group_by(Population) %>%
  summarise(
    mean_AF = mean(AF, na.rm = TRUE),
    sd_AF   = sd(AF,   na.rm = TRUE),
    n       = sum(!is.na(AF)),
    se_AF   = sd_AF / sqrt(n)
  )

ggplot(AF_summary, aes(x = Population, y = mean_AF)) +
  geom_col(width = 0.6, fill = "steelblue") +   # single colour for all bars
  geom_errorbar(
    aes(ymin = mean_AF - se_AF, ymax = mean_AF + se_AF),
    width = 0.15,
    color = "black"                             # error bars in black
  ) +
  labs(
    x = "Population",
    y = "Mean allele frequency"
  ) +
  coord_cartesian(ylim = c(0, max(AF_summary$mean_AF) * 1.3)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "none",
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(color = "black")
  )

#Export cleaned data into csv 
write.csv(AF_stats, "AF_stats_table.csv", row.names = FALSE)
