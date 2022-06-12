# Use the matched genes and merge back to the regions in each 
# DESeq2 output. 

import pandas as pd 

def main(matched_genes_file, deseq2_output): 
    matched = pd.read_csv(matched_genes_file, sep = '\t')
    
    return None 