# specify the input file path
input_file <- "../E.-coli-RBS-Prediction-and-Optimization/30nt_upstream_30nt_ORF.fasta"

# read in the multi-gene FASTA file
seqs <- readLines(input_file)

# split the sequences into individual sequences
seqs <- split(seqs, cumsum(grepl("^>", seqs)))

# iterate over each sequence and run RNAfold
for (i in seq_along(seqs)) {

  # extract the sequence name and sequence itself
  seq_name <- gsub(">", "", seqs[[i]][1])
  seq <- paste(seqs[[i]][-1], collapse = "")

  # create a temporary file to write the sequence
  seq_file <- tempfile()
  writeLines(seq, seq_file)

  # run RNAfold on the sequence file and capture the output
  fold_cmd <- paste0("RNAfold --noPS --noLP < ", seq_file)
  result <- system(fold_cmd, intern = TRUE)

  # extract the structure from the RNAfold output
  seq_structure <- result[1]

  # write the structure to a separate file for this sequence
  writeLines(seq_structure, paste0(seq_name, ".fold"))

  # delete the temporary sequence file
  file.remove(seq_file)
}
