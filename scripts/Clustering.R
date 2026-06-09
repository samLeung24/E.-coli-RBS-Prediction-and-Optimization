### Sam, misè a jour

library(Biostrings)
library(cluster)

# Define your RNA sequences
sequences <- readDNAStringSet("30nt_upstream_30nt_ORF.fasta")
# Get the number of sequences
n <- length(sequences)

# Create a pairwise distance matrix
distance_matrix <- matrix(0, nrow = n, ncol = n)

# Calculate pairwise distances using Needleman-Wunsch global alignment
for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      distance_matrix[i, j] <- 0
    } else {
      alignment <- pairwiseAlignment(sequences[i], sequences[j], type = "global", substitutionMatrix = nucleotideSubstitutionMatrix())
      distance_matrix[i, j] <- 1 - (alignment@score / max(width(sequences[i]), width(sequences[j])))
    }
  }
}

# Perform hierarchical clustering
hc <- hclust(as.dist(distance_matrix), method = "complete")

# Plot the dendrogram
pdf(file="hc.pdf")
plot(hc)
dev.off()
