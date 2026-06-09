table <- read.csv("forRF9.csv")
nRow <- length(table$Locus)
cai_tab <- data.frame(
  X6 <- substr(table$ORF,1,6),
  X12 <- substr(table$ORF,7,12),
  X18 <- substr(table$ORF,13,18),
  X24 <- substr(table$ORF,19,24),
  X30 <- substr(table$ORF,25,30)
)
names(cai_tab) <- c("X6","X12","X18","X24","X30")

write.csv(cai_tab,"cai_tab.csv",row.names = FALSE)

## After cai calculation with cai.py
cai_tab <- read.csv("cai_tab.csv")
cai_tab$avg_cai_ORF30 <- rowMeans(cai_tab[,c("score6","score12","score18",
                                          "score24","score30")])

table$avg_cai_ORF30 <- cai_tab$avg_cai_ORF30

write.csv(table,"forRF10.csv",row.names = FALSE)

gc_ORF30 <- substr(table$ORF,1,30)
table$gc_ORF30 <- str_count(gc_ORF30, "G|C")/nchar(gc_ORF30)
model <- lm(data = table,loged_el~gc_ORF30)
summary(model)

end <- substr(table$ORF,nchar(table$ORF)-49,nchar(table$ORF))
end_DNA <- DNAStringSet(end)
names(end_DNA) <- table$Locus
end_RNA <- RNAStringSet(end_DNA)

writeXStringSet(end_RNA,"../ORFEnd50.fasta")

structure <- c()
score <- c()

ORFEnd30_predict <- readLines("../ORFEnd50.predict")
for (i in 1:length(table$Locus)) {
  info <- strsplit(ORFEnd30_predict[3*i]," ")
  structure[i] <- info[[1]][1]
  score[i] <- gsub("[()]","",info[[1]][2])
}

table$ORFEnd50_structure <- structure
table$ORFEnd50_score <- as.numeric(score)

write.csv(table,"forRF11.csv",row.names = FALSE)
