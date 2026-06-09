library(stringr)
### Load positive genes from Prodigal SD identification
positives <- read.csv("bulkProdigal.csv")
GeneID <- positives[,"gene_id"]
Locus <- positives[,"locus_name"]
SDs <- positives[complete.cases(positives),"sd_seq"]
SDstarts <- positives[complete.cases(positives),"sd_start"]
SDends <- positives[complete.cases(positives),"sd_end"]

file1 <- readLines("../30_30.predict")
file2 <- readLines("../40_30.predict")
file3 <- readLines("../50_30.predict")

table <- data.frame(matrix(ncol = 11, nrow = length(GeneID)))
colnames(table) <- c("gene_id","sd_seq","30nt","30nt_proportion_paired", "30nt_score",
                     "40nt","40nt_proportion_paired", "40nt_score",
                     "50nt","50nt_proportion_paired", "50nt_score")

for (i in 1:length(GeneID)) {
  for (j in seq(from=1, to=length(file1),by=3)) {
    GI <- substring(sub("\\|.*", "",file1[j]),2)
    if (grepl(GeneID[i],GI,fixed = TRUE)) {
      table$gene_id[i] <- GeneID[i]
      startpose1 <- 0
      startpose2 <- 0
      startpose3 <- 0
      endpose1 <- 60
      endpose2 <- 70
      endpose3 <- 80

      table$`30nt`[i] <- substring(file1[j+2],first = startpose1, last = endpose1)
      table$`30nt_proportion_paired`[i] <- (sum(strsplit(table$`30nt`[i], "")[[1]] == "(")+
                                              sum(strsplit(table$`30nt`[i], "")[[1]] == ")"))/nchar(table$`30nt`[i])
      table$`30nt_score`[i] <- gsub("\\(|\\)", "",str_extract(file1[j+2],"\\(([^()]+)\\)[^()]*$"))
      table$`40nt`[i] <- substring(file2[j+2],first = startpose1, last = endpose1)
      table$`40nt_proportion_paired`[i] <- (sum(strsplit(table$`40nt`[i], "")[[1]] == "(")+
                                              sum(strsplit(table$`40nt`[i], "")[[1]] == ")"))/nchar(table$`30nt`[i])
      table$`40nt_score`[i] <- gsub("\\(|\\)", "",str_extract(file2[j+2],"\\(([^()]+)\\)[^()]*$"))
      table$`50nt`[i] <- substring(file3[j+2],first = startpose1, last = endpose1)
      table$`50nt_proportion_paired`[i] <- (sum(strsplit(table$`50nt`[i], "")[[1]] == "(")+
                                              sum(strsplit(table$`50nt`[i], "")[[1]] == ")"))/nchar(table$`30nt`[i])
      table$`50nt_score`[i] <- gsub("\\(|\\)", "",str_extract(file3[j+2],"\\(([^()]+)\\)[^()]*$"))
    }
  }
}

write.csv(table,file ="WholeProportionPairedMxfold2.csv",row.names = FALSE)
