library(stringr)

fasta_sequences <- read.csv("~/Desktop/projectRBS/prodigalPrediction.csv",sep=";",header = FALSE)
fasta_sequences <- fasta_sequences[,-14]
fasta_sequences$gene_id <- sub("\\|.*", "",fasta_sequences$V1)
fasta_sequences$locus_name <- sub("^[^|]+\\|([^,]+),.*", "\\1", fasta_sequences$V1)
fasta_sequences$sd_seq <- substring(fasta_sequences$V4, first = regexpr("=", fasta_sequences$V4) + 1)
fasta_sequences$spacer_length <- str_extract(fasta_sequences$V5, "(?<=\\=)\\d+|None")
fasta_sequences$gc_content <- substring(fasta_sequences$V6, first = regexpr("=", fasta_sequences$V6) + 1)

table <- read.csv("~/Desktop/projectRBS/E.-coli-RBS-Prediction-and-Optimization/predictionBioString4.csv")

tab_new <- fasta_sequences[,c("gene_id","sd_seq","spacer_length","gc_content")]
tab_need <- table[,c("gene_id","locus_name","RBS","ORF")]

merged_table <- merge(tab_new, tab_need, by = "gene_id", all.x = TRUE)

# initialize sd_start and sd_end columns with NA
merged_table$sd_start <- NA
merged_table$sd_end <- NA

# loop through each row
for (i in 1:nrow(merged_table)) {
  # check if sd_seq is not "None"
  if (merged_table$sd_seq[i] != "None") {
    # get start position of sd_seq within RBS
    start_pos <- regexpr(merged_table$sd_seq[i], merged_table$RBS[i])
    # check if sd_seq is present in RBS
    if (start_pos != -1) {
      # get end position of sd_seq within RBS
      end_pos <- start_pos + attr(regexpr(merged_table$sd_seq[i], merged_table$RBS[i]), "match.length") - 1
      # assign start and end positions to respective columns
      merged_table$sd_start[i] <- start_pos
      merged_table$sd_end[i] <- end_pos
    }
  }
}

merged_table$spacer_start <- NA
merged_table$spacer_end <- NA
merged_table$spacer_seq <- NA


for (i in seq(nrow(merged_table))) {
  # Get the RBS and sd_seq for this row
  rbs <- merged_table[i, "RBS"]
  sd_seq <- merged_table[i, "sd_seq"]
  
  # If sd_seq is "None", set spacer_start, spacer_end, and spacer_seq to NA
  if (sd_seq == "None") {
    merged_table[i, "spacer_start"] <- NA
    merged_table[i, "spacer_end"] <- NA
    merged_table[i, "spacer_seq"] <- NA
  } else {
    # Find the start and end positions of the sd_seq within the RBS
    sd_start <- regexpr(sd_seq, rbs)[1]
    sd_end <- sd_start + attr(regexpr(sd_seq, rbs), "match.length") - 1
    
    # Calculate the start and end positions of the spacer
    spacer_start <- sd_end + 1
    spacer_end <- nchar(rbs)
    
    # Extract the spacer sequence
    spacer_seq <- substr(rbs, spacer_start, spacer_end)
    
    # Update the merged table with the spacer start, end, and sequence
    merged_table[i, "spacer_start"] <- spacer_start
    merged_table[i, "spacer_end"] <- spacer_end
    merged_table[i, "spacer_seq"] <- spacer_seq
  }
}

final_table <- merged_table[, c("gene_id","locus_name","gc_content","sd_start",
                                 "sd_end", "sd_seq", "spacer_start",
                                 "spacer_end", "spacer_length", "spacer_seq",
                                 "RBS", "ORF")]

final_table[final_table == "None"] <- NA
percent_non_na <- sum(complete.cases(final_table)) / nrow(final_table) * 100
print(percent_non_na)

# Write the merged table to a CSV file
write.csv(final_table,"/Users/samleung/Desktop/projectRBS/E.-coli-RBS-Prediction-and-Optimization/bulkProdigal.csv", row.names = FALSE)
