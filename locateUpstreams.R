library(Biostrings)
library(seqinr)

genes <- readDNAStringSet("../NC_000913.3.fasta")
genome <- readDNAStringSet("../Sequence.fasta")
anti_genome <- complement(genome)
rev_genome <- reverse(anti_genome)
info <- read.csv("../processed_sequences_pax_dna_prot_Ecoli-K12_20230118.csv")
strand <- info$Strand

upstream_len <- 50
upstream_seqs <- DNAStringSet()
start_poses_bank <- integer()

for (i in 1:length(genes)) {
  
  # extract the 21st to 50th nucleotides of the gene sequence
  gene_seq <- substring(genes[i], 21, 50)
  
  # find all the location(s) of the gene in the genome, check existence
  gene_posf <-  vmatchPattern(gene_seq, genome)
  gene_posr <- vmatchPattern(gene_seq,rev_genome)
  start_poses <- integer()
  
  for (j in 1:length(start(gene_posf[[1]]))){
    start_poses <- c(start_poses,start(gene_posf[[1]])[j])
  }
  
  for (k in 1:length(start(gene_posr[[1]]))){
    temp_sp <- start(gene_posr[[1]])[k]
    real_sp <- width(genome)-temp_sp+1
    start_poses <- c(start_poses,real_sp)
  }
  
  # sort the start_poses
  start_poses <- sort(start_poses)
  
  # find the current start_pos to take
  start_pos <- start_poses[1]
  
  for (l in 1:100) {
   if (start_poses[l] %in% start_poses_bank & !i %in% c(458,473)) {
     start_pos <- start_poses[l+1]
   } else {
     start_poses_bank <- c(start_poses_bank,start_pos)
     break
   }
  }
  
  # check if start_pos is still empty
  if (is.na(start_pos)) {
    # skip this gene
    next
  }
  
  # calculate the position of upstream_start
  if (strand[i] == "+") {
    upstream_start <- start_pos - upstream_len
  } else {
    upstream_start <- start_pos + upstream_len
    pseudo_us <- width(genome)-upstream_start+1
  }

  # check if the upstream start position is less than 1
  if (!is.na(upstream_start)) {
    if (strand[i] == "+") {
      upstream_seq <- substring(genome, first=upstream_start,last=start_pos-1)
    } else {
      pseudo_start_pos <- width(genome)-start_pos+1
      upstream_seq <- substring(rev_genome, first=pseudo_us,last=pseudo_start_pos-1)
    }
  }
  
  # add the upstream sequence to the DNAStringSet
  upstream_seqs <- c(upstream_seqs, upstream_seq)
  
  # add the description to the upstream sequence
  names(upstream_seqs)[i] <- names(genes)[i]
}

tail(upstream_seqs)

write.fasta(as.list(upstream_seqs),names = names(upstream_seqs),file.out = "50nt_upstream.fasta")
