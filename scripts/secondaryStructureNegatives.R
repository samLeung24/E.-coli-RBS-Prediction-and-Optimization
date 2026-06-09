### Load negative genes from Prodigal SD identification
negatives <- read.csv("bulkProdigal.csv")
GeneID <- negatives[!complete.cases(positives),"gene_id"]
Locus <- negatives[!complete.cases(positives),"locus_name"]
SDs <- negatives[!complete.cases(positives),"sd_seq"]
SDstarts <- negatives[!complete.cases(positives),"sd_start"]
SDends <- negatives[!complete.cases(positives),"sd_end"]


dir1 <- "../30nt_2nd/"
dir2 <- "../40nt_2nd/"
dir3 <- "../50nt_2nd/"

tablen <- data.frame(matrix(ncol = 11, nrow = length(GeneID)))
colnames(tablen) <- c("gene_id","sd_seq","30nt","30nt_proportion_paired", "30nt_energy",
                     "40nt","40nt_proportion_paired", "40nt_energy",
                     "50nt","50nt_proportion_paired", "50nt_energy")



### Getting the base-pair status of SD in 2nd structure prediction
for (i in 1:length(GeneID)) {
  tablen$gene_id[i] <- GeneID[i]
  tablen$sd_seq[i] <- SDs[i]
  file1 <- readLines(paste0(dir1,GeneID[i],"|",Locus[i],".fold"), n = 2)
  file2 <- readLines(paste0(dir2,GeneID[i],"|",Locus[i],".fold"), n = 2)
  file3 <- readLines(paste0(dir3,GeneID[i],"|",Locus[i],".fold"), n = 2)

  tablen$`30nt_energy`[i] <- as.numeric(sub(".*\\((.*)\\)", "\\1", file1[2]))
  tablen$`40nt_energy`[i] <- as.numeric(sub(".*\\((.*)\\)", "\\1", file2[2]))
  tablen$`50nt_energy`[i] <- as.numeric(sub(".*\\((.*)\\)", "\\1", file3[2]))
}

t.test(table$`40nt_energy`,tablen$`40nt_energy`)

mean(table$`30nt_energy`)
mean(tablen$`30nt_energy`)
