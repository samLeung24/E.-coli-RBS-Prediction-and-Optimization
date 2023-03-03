### C'est Jeudi! One to go!

library(Biostrings)

# define the SD sequences
sd_template1 <- DNAString("AGG")
sd_template2 <- DNAString("GGA")
sd_template3 <- DNAString("GAG")

# read in the fasta file
fasta_sequences <- readDNAStringSet("~/Desktop/projectRBS/NC_000913.3.fasta")

# extract gene ID and locus name from fasta headers
header_fields <- strsplit(names(fasta_sequences), "\\|")
gene_ids <- sapply(header_fields, "[[", 1)
locus_names <- sapply(header_fields, "[[", 2)

# initialize variables for output
sd_start <- rep(NA, length(fasta_sequences))
sd_end <- rep(NA, length(fasta_sequences))
sd_seq <- rep(NA, length(fasta_sequences))
spacer_start <- rep(NA, length(fasta_sequences))
spacer_end <- rep(NA, length(fasta_sequences))
spacer_length <- rep(NA, length(fasta_sequences))
spacer_seq <- rep(NA, length(fasta_sequences))
RBS <- rep(NA, length(fasta_sequences))
ORF <- rep(NA, length(fasta_sequences))

# loop through each sequence in the fasta file
for (i in 1:length(fasta_sequences)) {
  seq <- fasta_sequences[i]
  seq_length <- length(seq)
  
  # adjustable matchmaking threshold
  sd_matches1 <- matchPattern(sd_template1, substring(seq, 1, 20), max.mismatch = 0, min.mismatch = 0)
  sd_matches2 <- matchPattern(sd_template2, substring(seq, 1, 20), max.mismatch = 0, min.mismatch = 0)
  sd_matches3 <- matchPattern(sd_template3, substring(seq, 1, 20), max.mismatch = 0, min.mismatch = 0)
  
  # if SD sequence is found, adjustable spacer threshold
  if (length(sd_matches1) > 0 & (20 - end(sd_matches1)[1] == 6 | 
                                20 - end(sd_matches1)[1] == 7 |
                                20 - end(sd_matches1)[1] == 8 |
                                20 - end(sd_matches1)[1] == 9)){
    # get start and end positions of SD sequence
    sd_start[i] <- start(sd_matches1)[1]
    sd_end[i] <- end(sd_matches1)[1]
    sd_seq[i] <- substring(seq,start(sd_matches1)[1],end(sd_matches1)[1])
    
    # define start and end positions of spacer and the spacer_seq
    spacer_start[i] <- sd_end[i] + 1
    spacer_end[i] <- 20
    spacer_seq[i] <- substring(seq,sd_end[i] + 1, 20)
    RBS[i] <- substring(seq, 1, 20)
    ORF[i] <- substring(seq,first = 21)
    
    # calculate spacer length
    spacer_length[i] <- spacer_end[i] - spacer_start[i] + 1
  } else if (length(sd_matches2) > 0 & (20 - end(sd_matches2)[1] == 6 | 
                                        20 - end(sd_matches2)[1] == 7 |
                                        20 - end(sd_matches2)[1] == 8 |
                                        20 - end(sd_matches2)[1] == 9)){
    # get start and end positions of SD sequence
    sd_start[i] <- start(sd_matches2)[1]
    sd_end[i] <- end(sd_matches2)[1]
    sd_seq[i] <- substring(seq,start(sd_matches2)[1],end(sd_matches2)[1])
    
    # define start and end positions of spacer and the spacer_seq
    spacer_start[i] <- sd_end[i] + 1
    spacer_end[i] <- 20
    spacer_seq[i] <- substring(seq,sd_end[i] + 1, 20)
    RBS[i] <- substring(seq, 1, 20)
    ORF[i] <- substring(seq,first = 21)
    
    # calculate spacer length
    spacer_length[i] <- spacer_end[i] - spacer_start[i] + 1
  } else if (length(sd_matches3) > 0 & (20 - end(sd_matches3)[1] == 6 | 
                                        20 - end(sd_matches3)[1] == 7 |
                                        20 - end(sd_matches3)[1] == 8 |
                                        20 - end(sd_matches3)[1] == 9)){
    sd_start[i] <- start(sd_matches3)[1]
    sd_end[i] <- end(sd_matches3)[1]
    sd_seq[i] <- substring(seq,start(sd_matches3)[1],end(sd_matches3)[1])
    
    # define start and end positions of spacer and the spacer_seq
    spacer_start[i] <- sd_end[i] + 1
    spacer_end[i] <- 20
    spacer_seq[i] <- substring(seq,sd_end[i] + 1, 20)
    RBS[i] <- substring(seq, 1, 20)
    ORF[i] <- substring(seq,first = 21)
    
    # calculate spacer length
    spacer_length[i] <- spacer_end[i] - spacer_start[i] + 1
  } else {
    # if SD sequence is not found, fill with NAs
    sd_start[i] <- NA
    sd_end[i] <- NA
    sd_seq[i] <- NA
    spacer_start[i] <- NA
    spacer_end[i] <- NA
    spacer_length[i] <- NA
    spacer_seq[i] <- NA
    RBS[i] <- substring(seq, 1, 20)
    ORF[i] <- substring(seq,first = 21)
  }
}

output_df <- data.frame(gene_id = gene_ids,
                        locus_name = locus_names,
                        sd_start = sd_start,
                        sd_end = sd_end,
                        sd_seq = sd_seq,
                        spacer_start = spacer_start,
                        spacer_end = spacer_end,
                        spacer_length = spacer_length,
                        spacer_seq = spacer_seq,
                        RBS = RBS,
                        ORF = ORF)

# write output to csv
write.csv(output_df, file = "predictionBioString4.csv", row.names = FALSE)

# test percent positive
percent_non_na <- sum(complete.cases(output_df)) / nrow(output_df) * 100
print(percent_non_na)
