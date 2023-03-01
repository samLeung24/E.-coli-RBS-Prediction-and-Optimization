### Comment tu t'applle aujourd'hui? C'est mecredi

library(Biostrings)

# define the SD sequence
sd_template <- DNAString("AGGA")
sd_template2 <- DNAString("AGGT")

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
spacer_seq <- rep(NA, length(fasta_sequences))
`5'UTR` <- rep(NA, length(fasta_sequences))
ORF <- rep(NA, length(fasta_sequences))

# loop through each sequence in the fasta file
for (i in 1:length(fasta_sequences)) {
  seq <- fasta_sequences[i]
  seq_length <- length(seq)
  
  # search for SD sequence in first 20 nucleotides
  sd_matches <- matchPattern(sd_template, substring(seq, 1, 20), max.mismatch = 1, min.mismatch = 0)
  sd_matches2 <- matchPattern(sd_template2, substring(seq, 1, 20), max.mismatch = 1, min.mismatch = 0)
  # if SD sequence is found
  if (length(sd_matches) > 0) {
    # get start and end positions of SD sequence
    sd_start[i] <- start(sd_matches)[1]
    sd_end[i] <- end(sd_matches)[1]
    
    # define start and end positions of spacer and the spacer_seq
    spacer_start[i] <- sd_end[i] + 1
    spacer_end[i] <- 20
    spacer_seq[i] <- substring(seq,sd_end[i] + 1, 20)
    `5'UTR`[i] <- substring(seq, 1, 20)
    ORF[i] <- substring(seq,first = 21)
    
    # calculate spacer length
    spacer_length[i] <- spacer_end[i] - spacer_start[i] + 1
  } else if (length(sd_matches2) > 0){
    # if initial matching does not yield
    # get start and end positions of SD sequence
    sd_start[i] <- start(sd_matches2)[1]
    sd_end[i] <- end(sd_matches2)[1]
    
    # define start and end positions of spacer and the spacer_seq
    spacer_start[i] <- sd_end[i] + 1
    spacer_end[i] <- 20
    spacer_seq[i] <- substring(seq,sd_end[i] + 1, 20)
    `5'UTR`[i] <- substring(seq, 1, 20)
    ORF[i] <- substring(seq,first = 21)
    
  } else {
    # if SD sequence is not found, fill with NAs
    sd_start[i] <- NA
    sd_end[i] <- NA
    spacer_start[i] <- NA
    spacer_end[i] <- NA
    spacer_length[i] <- NA
    spacer_seq[i] <- NA
    `5'UTR`[i] <- substring(seq, 1, 20)
    ORF[i] <- substring(seq,first = 21)
  }
}

# create data frame for output
output_df <- data.frame(gene_id = gene_ids,
                        locus_name = locus_names,
                        sd_start = sd_start,
                        sd_end = sd_end,
                        spacer_start = spacer_start,
                        spacer_end = spacer_end,
                        spacer_length = spacer_length,
                        spacer_seq = spacer_seq,
                        `5'UTR` = `5'UTR`,
                        ORF = ORF)

# write output to csv
write.csv(output_df, file = "predictionBioString1.csv", row.names = FALSE)
