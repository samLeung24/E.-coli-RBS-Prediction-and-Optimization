library(dplyr)

table <- read.csv("forRF8.csv")

output <- c()
for (i in 1:length(table$ORF)){
  ncodon_each <- nchar(table$ORF[i])/3
  for (j in 1:ncodon_each) {
    output[[i]][j]<- substr(table$ORF[i],3*j-2,3*j)
  }
}
output$Locus <- table$Locus
Codon_Genes <- output %>%
  relocate(Locus)


Codon <- c("GCT","GCC","GCG","GCA",
                 "GGT","GGC","GGG","GGA",
                 "CCT","CCC","CCG","CCA",
                 "ACT","ACC","ACG","ACA",
                 "GTT","GTC","GTG","GTA",
                 "TCT","TCC","TCG","TCA","AGT","AGC",
                 "CGT","CGC","CGG","CGA","AGG","AGA",
                 "CTT","CTC","CTG","CTA","TTG","TTA",
                 "TTT","TTC",
                 "AAT","AAC",
                 "AAG","AAA",
                 "GAT","GAC",
                 "GAG","GAA",
                 "CAT","CAC",
                 "CAG","CAA",
                 "ATT","ATC","ATA",
                 "ATG",
                 "TAT","TAC",
                 "TGA","TAG","TAA",
                 "TGT","TGC",
                 "TGG")
Percent <- as.numeric(c(1.52,2.55,3.37,2.02,2.47,2.96,1.11,0.79,0.7,0.55,2.33,0.85,
             0.89,2.34,1.44,0.71,1.83,1.53,2.62,1.09,0.8,0.86,0.89,0.72,0.87,1.6,
             2.09,2.2,0.54,0.36,0.12,0.21,1.1,1.11,5.28,0.39,1.37,1.39,
             2.22,1.65,1.77,2.15,1.03,3.37,3.21,1.91,1.78,3.95,1.29,0.97,2.89,1.54,
             3.04,2.52,0.44,2.78,1.6,1.22,0.11,0.03,0.21,0.5,0.65,1.53))

Codon_Usage <- as.data.frame(cbind(Codon,Percent))
Codon_Usage$Percent <- as.numeric(Codon_Usage$Percent)/100

Codon_Score <- data.frame(matrix(0, nrow = 3911, ncol = 2359))
for (i in 1:length(Codon_Genes$Locus)) {
  for (j in 2:length(Codon_Genes)) {
    if (is.na(Codon_Genes[i,j])) {
      for(k in j:length(Codon_Genes)) {
        Codon_Score[i,j] <- NA
      }
      break
    }
    Codon_Score[i,j] <- Codon_Usage[Codon == Codon_Genes[i,j],2]
  }
}
Codon_Score[,1] <- Codon_Genes$Locus
names(Codon_Score)[1] <- "Locus"

Codon_Preference <- data.frame(Locus = Codon_Score$Locus,
                               Avg_Score = vector("numeric", length(Codon_Score$Locus)))
for (i in 1:length(Codon_Preference$Locus)) {
  stop_pos = 2360
  for (j in 2:length(Codon_Score)) {
    if (is.na(Codon_Score[i,j])) {
      stop_pos = j-1
      break
    }
  }
  Codon_Preference$Avg_Score[i] <- mean(unlist(Codon_Score[i,2:stop_pos]))
}

new_table <- merge(table,Codon_Preference,by = "Locus")
new_table <- distinct(new_table, Locus, .keep_all = TRUE)


### Alternative method (now in used)
table <- read.csv("forRF8.csv")
Codon_Genes <- vector("list", length = length(table$Locus))
names(Codon_Genes) <- table$Locus
for (i in 1:length(Codon_Genes)) {
  ncodon_each <- nchar(table$ORF[i])/3
  for (j in 1:ncodon_each) {
    Codon_Genes[[i]][j] <- substr(table$ORF[i],3*j-2,3*j)
  }
}

### Info of all tRNA and corresponding
Codon_ref2 <- read.table("../all_tRNA_info.txt",skip = 3, sep = "\t")

All_genes <- read.csv("../All-genes-of-E.-coli-K-12-substr.csv")
All_genes$Locus <- substr(All_genes$Sequence...coordinates.of.DNA.region,2,5)
Ref_finale <- data.frame("Anticodon" = vector("character",length = 87),
                    "Locus" = vector("character",length = 87))

for (i in 1:length(Codon_ref2$V1)) {
  Ref_finale$Anticodon[i] <- Codon_ref2$V6[i]
  if (Codon_ref2$V3[i] %in% All_genes$Left.End.Position) {
    Ref_finale$Locus[i] = All_genes[which(All_genes$Left.End.Position == Codon_ref2$V3[i]), "Locus"]
  } else if (Codon_ref2$V3[i] %in% All_genes$Right.End.Position) {
    Ref_finale$Locus[i] = All_genes[which(All_genes$Right.End.Position == Codon_ref2$V3[i]), "Locus"]
  } else if (Codon_ref2$V4[i] %in% All_genes$Left.End.Position) {
    Ref_finale$Locus[i] = All_genes[which(All_genes$Left.End.Position == Codon_ref2$V4[i]), "Locus"]
  } else if (Codon_ref2$V4[i] %in% All_genes$Right.End.Position) {
    Ref_finale$Locus[i] = All_genes[which(All_genes$Right.End.Position == Codon_ref2$V4[i]), "Locus"]
  } else Ref_finale$Locus[i] = NA
}

Codon_Available <- read.csv("~/Desktop/projectRBS/Supplemental File S6.csv")
Codon_Available$mean <- rowMeans(Codon_Available[2:10])
sum <- sum(Codon_Available$mean)
Codon_Available$`#/1000` <- 1000*(Codon_Available$mean/sum)
Codon_Available <- merge.data.frame(Ref_finale,Codon_Available, by.x = "Locus", by.y = "tRNA", all = TRUE)
Codon_Available[35,"Anticodon"] <- "GTG"
Codon_Available <- Codon_Available[complete.cases(Codon_Available),]

base_compliment <- function(base) {
  base_c <- "N"
  if (base == "A") base_c <- "T"
    else if (base == "T") base_c <- "A"
    else if (base == "C") base_c <- "G"
    else if (base == "G") base_c <- "C"
    return(base_c)
}

for (i in 1:length(Codon_Available$Anticodon)) {
  Codoni <- Codon_Available$Anticodon[i]
  base1 <- substr(Codoni,1,1)
  base2 <- substr(Codoni,2,2)
  base3 <- substr(Codoni,3,3)
  Codon_Available$Codon[i] <- paste0(base_compliment(base3),
                                     base_compliment(base2),
                                     base_compliment(base1))
}

Codon_Score2 <- vector("list", length = length(table$Locus))
names(Codon_Score2) <- table$Locus
for (i in 1:length(Codon_Genes)) {
  for (j in 1:length(Codon_Genes[[i]])) {
    print(c(i,j))
    if (Codon_Genes[[i]][j] %in% Codon_Available$Codon)
    Codon_Score2[[i]][j] <- Codon_Available[Codon_Available$Codon == Codon_Genes[[i]][j],"#/1000"]
    else Codon_Score2[[i]][j] <- 0
  }
}

Codon_Preference2 <- data.frame("Codon" = names(Codon_Score2),
                                "Avg_Score2" = vector("character",length = length(Codon_Score2)))
for (i in 1:length(Codon_Genes)) {
  print(i)
  Codon_Preference2$Avg_Score2[i] <- mean(Codon_Score2[[i]])
}

table$Avg_Score2 <- as.numeric(Codon_Preference2$Avg_Score2)

write.csv(table,"forRF8.csv",row.names = FALSE)
