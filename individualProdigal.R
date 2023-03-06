library(tidyverse)
library(dplyr)

# Get all the individual CSV files
csv_files <- list.files(pattern = "*.csv", full.names = FALSE)

# Read in the first CSV file to create the merged table
mod_table <- read.csv(csv_files[2], sep = ";", header = FALSE)
vector <- sub("=.*", "", mod_table[1,])
merged_table <- data.frame(matrix(ncol = length(vector)+2))
colnames(merged_table) <- c("gene_id","locus_name",vector)

# Loop through the remaining CSV files and merge them with the first one
for (i in 1:length(csv_files)) {
  file_name <- csv_files[i]
  file_parts <- strsplit(file_name, "|", fixed = TRUE)[[1]]
  gene_id <- file_parts[1]
  locus_name <- file_parts[2]
  
  # If the individual table has no rows, fill NAs and continue to the next file
  if (file.info(csv_files[i])$size == 0) {
    individual_table <- data.frame(matrix(NA, nrow = 1, ncol = ncol(merged_table)))
    colnames(individual_table) <- c("gene_id","locus_name",vector)
    individual_table$gene_id <- gene_id
    individual_table$locus_name <- locus_name
  } else {
    temp_table <- data.frame(matrix(ncol = 16))
    individual_table <- read.csv(csv_files[i], sep = ";", header = FALSE)
    data <- sub(".*=", "", individual_table[1,])
    colnames(temp_table) <- c("gene_id","locus_name",vector)
    temp_table[1,] <- c(gene_id,locus_name,data)
    individual_table <- temp_table
  }
  merged_table <- rbind(merged_table, individual_table)
}

final_table <- merged_table[-1,-16]
final_table <- final_table %>% 
  arrange(as.numeric(gene_id))

final_table[final_table == "None"] <- NA
percent_non_na <- sum(complete.cases(final_table)) / nrow(final_table) * 100
print(percent_non_na)

# Write the merged table to a CSV file
write.csv(final_table,"/Users/samleung/Desktop/projectRBS/E.-coli-RBS-Prediction-and-Optimization/individualProdigal.csv", row.names = FALSE)
