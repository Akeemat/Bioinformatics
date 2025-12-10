### ClinVar
clinvar_data <- read.delim("clinvar_result.txt")

#confirm file details
head(clinvar_data)
str(clinvar_data)
nrow(clinvar_data)
table(clinvar_data$Condition.s.)

# Deduplicate and filter for relevant data
clinvar_data_unique <- clinvar_data %>%
  distinct() %>%                                    # remove exact duplicates
  filter(str_detect(Gene.s., "ACE2")) %>%           # keep only ACE2 rows
  filter(!str_detect(str_to_lower(Condition.s.),
                     "not provided|not specified|see cases"))            

nrow(clinvar_data_unique)

#Generate summary table for clinical conditions only
clinvar_data_summary <- as.data.frame(table(clinvar_data_unique$Condition.s.))
colnames(clinvar_data_summary) <- c("Clinical condition", "Frequency")
clinvar_data_summary <- clinvar_data_summary[order(clinvar_data_summary$Frequency, decreasing = TRUE), ]
sum(clinvar_data_summary$Frequency)

# Generate a Cross-tab of variant types vs clinical conditions
variant_condition_summary <- clinvar_data_unique %>%
  group_by(Variant.type, Condition.s.) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  arrange(desc(Frequency))
view(variant_condition_summary)

#Export cleaned data into csv 
write.csv(clinvar_data_unique, "clinvar_data_unique.csv", row.names = FALSE)