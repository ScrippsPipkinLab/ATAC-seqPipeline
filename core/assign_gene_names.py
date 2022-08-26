import sys 
import pandas as pd
import os
# Take in merged BED file and assign gene names to each row from reference gtf 

# matched_genes_file = sys.argv[1]
# deseq2_folder = sys.argv[2]
# output_file = sys.argv[3]

matched_genes_file = '/blue/m.pipkin/stsuda/ATAC_Ets1KO_invitro/merged/matched_genes.txt'
deseq2_folder = '/blue/m.pipkin/stsuda/ATAC_Ets1KO_invitro/deseq2'
output_folder = '/blue/m.pipkin/stsuda/ATAC_Ets1KO_invitro/gene_labeled'


def closest_annotation(matched_genes_file, deseq2_output, output_file):
    '''
    Use Bedtools closest match to find closest annotation to each merged BED row
    '''
    matched_headers = [
        'Chr', 'Start', 'End', 'Strand', 'Score', 
        'Chr_closest', 'Start_closest', 'End_closest', 
        'Gene_closest', 'empty', 'Strand_closest',
        'Source_closest', 'type', 'empty_2', 
        'info_string', 'distance_closest'
        ]
    matched = pd.read_csv(matched_genes_file, sep='\t', header=None, names=matched_headers)
    matched = matched.astype({
        'Chr': str, 'Start': int, 'End': int, 'Strand': str, 'Score': float, 
        'Chr_closest': str, 'Start_closest': int, 'End_closest': int, 
        'Gene_closest': str, 'empty': str, 'Strand_closest': str,
        'Source_closest': str, 'type': str, 'empty_2': str, 
        'info_string': str, 'distance_closest': int
        })
        
    deseq_names = ['peak_id', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    deseq = pd.read_csv(deseq2_output, names=deseq_names, header=None, skiprows=1)
    # deseq = pd.read_csv(deseq2_output)
    print(deseq.head())
    split_peak_ids = deseq['peak_id'].str.split('_', n=4, expand=True)
    deseq['Chr'] = split_peak_ids[0] + '_' + split_peak_ids[1]
    deseq['Start'] = split_peak_ids[2]
    deseq['End'] = split_peak_ids[3]
    deseq = deseq.astype({
        'peak_id': str, 'baseMean': float, 'log2FoldChange': float, 
        'lfcSE': float, 'pvalue': float, 'padj': float, 
        'Chr': str, 'Start': int, 'End': int
        })

    merged = pd.merge(deseq, matched, how='left', on=['Chr', 'Start', 'End'])
    merged.to_csv(output_file, header=True)  

    return None 

# 
if __name__ == '__main__':
    for deseq2_output in os.listdir(deseq2_folder):
        if deseq2_output.endswith('.csv'):
            print(deseq2_output)
            print(deseq2_folder + '/' + deseq2_output)
            print(output_folder + '/' + deseq2_output)

            closest_annotation(
                matched_genes_file, 
                deseq2_folder + '/' + deseq2_output, 
                output_folder + '/' + deseq2_output
                )