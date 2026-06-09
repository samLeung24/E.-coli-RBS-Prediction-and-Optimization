library(Biostrings)

genes <- read.csv("forRF7.csv")
readRBS <- DNAStringSet(genes$RBS)
  RBS <- RNAStringSet(readRBS)
Genome <- readDNAStringSet("../MG1655_complete_genome.fasta")
  rRNA_start <- 4166659
  rRNA_end <- 4168200
  rRNA <- Genome[[1]][rRNA_start:rRNA_end]
  rRNA <- RNAStringSet(rRNA)

calculate_complementarity <- function(seq1,seq2,start,end) {
  ## Check for input validity
  if (!is(seq1, "RNAString") || !is(seq2, "RNAString")) {
    stop("Both input sequences must be RNAString objects.")
  }
  if (!is(start, "numeric") || !is(end, "numeric")) {
    stop("Both position indices must be numeric")
  }

  seq1 <- unlist(strsplit(as.character(seq1[start:end]),""))
  seq2 <- seq2[start:end]
  seq2_rev <-  reverse(rRNA)[[1]]
  seq2_comp <- unlist(strsplit(as.character(complement(seq2_rev)[1:20]),""))
  matches <- seq1 == seq2_comp
  score <- sum(matches)/(end+1-start)
  return(score)
}

scores <- vector(mode = "numeric", length(genes[,1]))

for (i in 1:length(genes[,1])) {
   score <- calculate_complementarity(RBS[[i]],rRNA[[1]],1,20)
   scores[i] <- score
}

genes$Complim <- scores

write.csv(genes,row.names = FALSE,file = "forRF8.csv")
