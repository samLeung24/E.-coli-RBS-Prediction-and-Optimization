### Load positive genes from Prodigal SD identification
positives <- read.csv("bulkProdigal.csv")
GeneID <- positives[complete.cases(positives),"gene_id"]
Locus <- positives[complete.cases(positives),"locus_name"]
SDs <- positives[complete.cases(positives),"sd_seq"]
SDstarts <- positives[complete.cases(positives),"sd_start"]
SDends <- positives[complete.cases(positives),"sd_end"]


dir1 <- "../30nt_2nd/"
dir2 <- "../40nt_2nd/"
dir3 <- "../50nt_2nd/"

table <- data.frame(matrix(ncol = 8, nrow = length(GeneID)))
colnames(table) <- c("gene_id","sd_seq","30nt","30nt_proportion_paired",
                     "40nt","40nt_proportion_paired","50nt","50nt_proportion_paired")


### Getting the base-pair status of SD in 2nd structure prediction
for (i in 1:length(GeneID)) {
  table$gene_id[i] <- GeneID[i]
  table$sd_seq[i] <- SDs[i]
  file1 <- readLines(paste0(dir1,GeneID[i],"|",Locus[i],".fold"), n = 2)
  file2 <- readLines(paste0(dir2,GeneID[i],"|",Locus[i],".fold"), n = 2)
  file3 <- readLines(paste0(dir3,GeneID[i],"|",Locus[i],".fold"), n = 2)
  startpose1 <- 10+SDstarts[i]
  startpose2 <- 20+SDstarts[i]
  startpose3 <- 30+SDstarts[i]
  endpose1 <- 10+SDends[i]
  endpose2 <- 20+SDends[i]
  endpose3 <- 30+SDends[i]
  
  table$`30nt`[i] <- substring(file1[2],first = startpose1, last = endpose1)
  table$`30nt_proportion_paired`[i] <- (sum(strsplit(table$`30nt`[i], "")[[1]] == "(")+
                                          sum(strsplit(table$`30nt`[i], "")[[1]] == ")"))/nchar(table$`30nt`[i])
  table$`40nt`[i] <- substring(file2[2],first = startpose2, last = endpose2)
  table$`40nt_proportion_paired`[i] <- (sum(strsplit(table$`40nt`[i], "")[[1]] == "(")+
                                          sum(strsplit(table$`40nt`[i], "")[[1]] == ")"))/nchar(table$`40nt`[i])
  table$`50nt`[i] <- substring(file3[2],first = startpose3, last = endpose3)
  table$`50nt_proportion_paired`[i] <- (sum(strsplit(table$`50nt`[i], "")[[1]] == "(")+
                                          sum(strsplit(table$`50nt`[i], "")[[1]] == ")"))/nchar(table$`50nt`[i])
}

paste("The mean proportion of paired SD nucleotides in 30nt-upstreams model is",round(mean(table$`30nt_proportion_paired`),3))
paste("The mean proportion of paired SD nucleotides in 40nt-upstreams model is",round(mean(table$`40nt_proportion_paired`),3))
paste("The mean proportion of paired SD nucleotides in 50nt-upstreams model is",round(mean(table$`50nt_proportion_paired`),3))


write.csv(table,"proportionPaired.csv")




