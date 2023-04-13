library(Biostrings)
library(seqinr)

upstream <- readDNAStringSet("40nt_upstream.fasta")
upAndORF <- readDNAStringSet("../NC_000913.3.fasta")
cut <- DNAStringSet()
sample <- DNAStringSet()

for (i in 1:length(upstream)) {
  cut[i] <- substring(upAndORF[i],21,50)
}

for (j in 1:length(upstream)) {
  sample[j] <- DNAStringSet(paste0(upstream[j],cut[j]))
  names(sample)[j] <- names(upstream)[j]
}


write.fasta(as.list(sample),names = names(sample),file.out = "40nt_upstream_30nt_ORF.fasta")

