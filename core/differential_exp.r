#!/usr/bin/env Rscript


library(DESeq2, quietly=TRUE)
library(stringr)
library(ashr)
library(EnhancedVolcano)
library(scales)

args = commandArgs(trailingOnly=TRUE)
cts = read.table(args[1], header=TRUE, row.names=1)
ssheet = read.table(args[2], header=TRUE, stringsAsFactors=TRUE)
results_dir = args[3]
# Assumes the results_dir has already been created by python 

####
coldata <- ssheet[, c("Status", "CT")]
rownames(coldata) <- ssheet[, "SampleName"]
print(coldata)


name.vec <- c("Chr", "Start", "End", "Strand", "Length")
for (i in 6:ncol(cts)) {
    name.vec <- append(name.vec, as.character(ssheet$SampleName[i-5]))
}
names(cts) <- name.vec


mat = cts[,6:ncol(cts)]
mat[is.na(mat)] <- 0

print(str(mat))
print(head(mat))

print(str(coldata))
print(head(coldata))

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = coldata,
                              design= ~ Status)
dds <- DESeq(dds)
resultsNames(dds)
#########
for (control in unique(ssheet$Status)) {
    for (treatment in unique(ssheet$Status)) {
        if (control == treatment){
            next
        } else {
            print(paste("Contrast: ", control, " vs ", treatment))
            res <- results(dds, contrast=c("Status", control, treatment),
                        independentFiltering=TRUE, parallel=TRUE)
            res <- lfcShrink(dds, contrast=c("Status", control, treatment), res=res, type="ashr")
            print(paste0(results_dir, "/", control, "_vs_", treatment, ".csv"))
            write.csv(as.data.frame(res), paste0(results_dir, "/", control, "_vs_", treatment, ".csv"))

            print(paste0(results_dir, "/", control, "_vs_", treatment, ".pdf"))
            print(paste(control, "vs.", treatment))
            volcano <- EnhancedVolcano(res,
                            lab = rownames(res),
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            title = paste(control, "vs.", treatment),
                            pCutoff = 0.05,
                            FCcutoff = 1.2
                            )
            ggsave(paste0(results_dir, "/", control, "_vs_", treatment, ".pdf"), volcano)

            # # Make pairwse peak/count plots for each sample
            control_samples <- as.character(ssheet[ssheet$Status == control,"SampleName"])
            treatment_samples <- as.character(ssheet[ssheet$Status == treatment,"SampleName"])
            # rowMeans behaes extremely stupidly and doesn't work on a dataframe with only one column 
            # This is why R sucks...
            if (length(control_samples) > 1) {
                control_counts <- rowMeans(mat[, control_samples])
            } else {
                control_counts <- mat[, control_samples]
            }
            if (length(treatment_samples) > 1) {
                treatment_counts <- rowMeans(mat[, treatment_samples])
            } else {
                treatment_counts <- mat[, treatment_samples]
            }
            ppc_df <- as.data.frame(cbind(control_counts, treatment_counts))
            ppc <- ggplot(ppc_df, aes(x=control_counts,
                                    y=treatment_counts),
                        ) +
                geom_density_2d() +
                geom_abline(intercept=0, slope=1, color="Black") +
                scale_y_continuous(trans='log10', 
                                    breaks = trans_breaks('log10', function(x) 10^x),
                                    labels = trans_format('log10', math_format(10^.x))
                                    ) +
                scale_x_continuous(trans='log10', 
                                    breaks = trans_breaks('log10', function(x) 10^x),
                                    labels = trans_format('log10', math_format(10^.x))
                                    ) +
                theme_minimal() +
                xlab(paste0(control)) +
                ylab(paste0(treatment)) +
                ggtitle("ATAC-seq mean counts/peak")
            ###### These plots may also be used
            # ppc <- ggplot(ppc_df, aes(x=control_counts,
            #                           y=treatment_counts),
            #               ) +
            #        geom_hex() +
            #        geom_abline(intercept=0, slope=1, color="Black") +
            #        scale_fill_continuous(type = "viridis") +
            #        scale_y_continuous(trans='log10') +
            #        scale_x_continuous(trans='log10')
            # ppc <- ggplot(ppc_df, aes(x=control_counts,
            #                           y=treatment_counts),
            #               ) +
            #        geom_point(alpha=0.05) +
            #        geom_abline(intercept=0, slope=1, color="Black") +
            #        scale_y_continuous(trans='log10') +
            #        scale_x_continuous(trans='log10')
            #######
            ggsave(paste0(results_dir, "/", control, "_vs_", treatment, "_ppc",".pdf"), ppc)
        }
    }
}

#########
main <- function(){
}
