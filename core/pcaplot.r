library(ggfortify)
library(stringr)
library(factoextra)

args = commandArgs(trailingOnly=TRUE)
cts = read.table(args[1], header=TRUE, row.names=1)
ssheet = read.table(args[2], header=TRUE, stringsAsFactors=TRUE)
results_dir = args[3]


print(str(cts))

counts <- cts[,6:ncol(cts)]
counts[is.na(counts)] <- 0


labels <- str_match(names(counts), "[a-zA-Z0-9-_]+\\.bam$")
print(labels)
names(counts) <- sub(".bam", "", labels)

###
counts$KO11 <- NULL


counts <- t(counts)
merged_df <- ssheet[order(match(ssheet, rownames(counts)))]
print(merged_df)

head(merged_df)
pca_prcomp <- prcomp(counts, center=TRUE, scale=TRUE)
print(summary(pca_prcomp))

pdf(paste0(results_dir, "/", "pca_elbow.pdf"))
screeplot(pca_prcomp, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
dev.off()


pdf(paste0(results_dir, "/", "pca_plot.pdf"))
biplot(pca_prcomp)
# autoplot(
#     pca_prcomp, label=TRUE, repel=TRUE, 
#     #colour=merged_df$Time, shape=merged_df$Treatment,
#     loadings.label=TRUE, loadings=TRUE
#     )
#legend("right")
dev.off()