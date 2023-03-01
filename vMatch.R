library(Biostrings)

# define the SD sequence
sd_template <- DNAString("AGGA")

# read in the fasta file
fasta_sequences <- readDNAStringSet("~/Desktop/projectRBS/NC_000913.3.fasta")

# extract gene ID and locus name from fasta headers
header_fields <- strsplit(names(fasta_sequences), "\\|")
gene_ids <- sapply(header_fields, "[[", 1)
locus_names <- sapply(header_fields, "[[", 2)

# initialize variables for output
sd_start <- rep(NA, length(fasta_sequences))
sd_end <- rep(NA, length(fasta_sequences))
spacer_start <- rep(NA, length(fasta_sequences))
spacer_end <- rep(NA, length(fasta_sequences))
spacer_length <- rep(NA, length(fasta_sequences))

# loop through each sequence in the fasta file
for (i in 1:length(fasta_sequences)) {
  seq <- fasta_sequences[i]
  seq_length <- length(seq)
  
  # search for SD sequence in first 20 nucleotides
  sd_matches <- matchPattern(sd_template, substring(seq, 1, 20), max.mismatch = 1, min.mismatch = 0)
  
  # if SD sequence is found
  if (length(sd_matches) > 0) {
    # get start and end positions of SD sequence
    sd_start[i] <- start(sd_matches)[1]
    sd_end[i] <- end(sd_matches)[1]
    
    # define start and end positions of spacer
    spacer_start[i] <- sd_end[i] + 1
    spacer_end[i] <- 20
    
    # calculate spacer length
    spacer_length[i] <- spacer_end[i] - spacer_start[i] + 1
  } else {
    # if SD sequence is not found, fill with NAs
    sd_start[i] <- NA
    sd_end[i] <- NA
    spacer_start[i] <- NA
    spacer_end[i] <- NA
    spacer_length[i] <- NA
  }
}

# create data frame for output
output_df <- data.frame(gene_id = gene_ids,
                        locus_name = locus_names,
                        sd_start = sd_start,
                        sd_end = sd_end,
                        spacer_start = spacer_start,
                        spacer_end = spacer_end,
                        spacer_length = spacer_length)

# write output to csv
write.csv(output_df, file = "NC_000913.3.csv", row.names = FALSE)
