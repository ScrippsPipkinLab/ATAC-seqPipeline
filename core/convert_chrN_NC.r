##########
# Convert Blacklist regions from mm10 to mm39
# Blacklist regions are from ENCODE 
# https://github.com/Boyle-Lab/Blacklist/tree/master/lists
# First liftOver using UCSC liftOver to mm39 
# Then convert chrN to pesky NC_ notation for Refseq 


lifted <- read.delim("/gpfs/home/snagaraja/ATACseqPipeline/refs/ENCODE_blacklist_mm10_liftover_m39.bed", header = FALSE)
names(lifted) <- c("chr", "start", "end")

report <- read.delim("/gpfs/home/snagaraja/ATACseqPipeline/refs/mm39_BuildReport.txt", comment.char = "#")
report$chrName <- paste("chr", report$Assigned.Molecule, sep="")

merged <- merge(lifted, report, by.x = "chr", by.y = "chrName")
final <- merged[, c("chr", "start", "end")]

write.table(final, "/gpfs/home/snagaraja/ATACseqPipeline/refs/blacklist.bed",
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

