import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import pandas as pd 
import argparse, sys, os

#### Command Line Args ####
parser = argparse.ArgumentParser(description='Draw venn diagrams from two lists of peaks')
parser.add_argument('-A', help='narrowPeak file for sample A', required=True)
parser.add_argument('-B', help='narrowPeak file for sample B', required=True)
parser.add_argument('-O', help='output directory path', required=True)
args = parser.parse_args()
#### ####


# SSHEET_PATH = '/blue/m.pipkin/stsuda/ATAC_Ets1KO_invitro/Ets1_ATAC.txt'

# # Read sample sheet file 
# ssheet = pd.read_csv(SSHEET_PATH, sep='\t')
# print(ssheet.head())
# # pseudocode: 
# # 1. convert bedgraph to list of peaks in the format [chr_start_end, chr_start_end, ...]
# # 2. Draw venn diagrams

def narrowPeak_to_list(narrowPeak_file):
    '''
    Convert narrowPeak file to list of peaks in the format [chr_start_end, chr_start_end, ...]
    '''
    sample_name = os.path.basename(narrowPeak_file)
    sample_name = sample_name.replace('_peaks.narrowPeak', '')
    peaks = []
    with open(narrowPeak_file) as f:
        for line in f:
            line = line.strip().split()
            peaks.append(line[0] + '_' + line[1] + '_' + line[2])
    
    return peaks, sample_name

# def draw_venns(peaks_list_A, peaks_list_B):
#     '''
#     Draw venn diagrams from two lists of peaks
#     '''

#     A_only = set(peaks_list_A) - set(peaks_list_B)
#     intersection = set(peaks_list_A) & set(peaks_list_B)
#     B_only = set(peaks_list_B) - set(peaks_list_A)


#     venn2(subsets=(len(A_only), len(intersection), len(B_only)), set_labels=('A', 'B'), alpha=0.5)
#     venn2_circles(subsets=(len(A_only), len(intersection), len(B_only)), linewidth=1.0)

#     plt.show()
#     plt.savefig('test.pdf')

#     return None

def main():

    # Set sample names and file paths to narrowpeaks 
    converted_peaks_A = narrowPeak_to_list(args.A)
    converted_peaks_B = narrowPeak_to_list(args.B)

    sample_A_name, sample_B_name = converted_peaks_A[1], converted_peaks_B[1]
    peaks_list_A, peaks_list_B = converted_peaks_A[0], converted_peaks_B[0]
    # Compute intersection and difference between A and B
    A_only = set(peaks_list_A) - set(peaks_list_B)
    intersection = set(peaks_list_A) & set(peaks_list_B)
    B_only = set(peaks_list_B) - set(peaks_list_A)
    total = len(A_only) + len(B_only) + len(intersection)

    print(f'{sample_A_name} only: {len(A_only)}')
    print(f'{sample_B_name} only: {len(B_only)}')
    print(f'{sample_A_name} and {sample_B_name} intersection: {len(intersection)}')
    # Draw venn diagrams
    venn2(subsets=(len(A_only), len(intersection), len(B_only)), set_labels=(sample_A_name, sample_B_name), alpha=0.5)
    venn2_circles(subsets=(len(A_only), len(intersection), len(B_only)), linewidth=1.0)
    print(f'Saving to: {args.O}/{sample_A_name}_and_{sample_B_name}.svg')
    # Save venn diagram to file
    plt.savefig(f'{args.O}/{sample_A_name}_and_{sample_B_name}.svg')
    plt.close()
    # Plot with Percentage labels 
    venn2(subsets=(len(A_only), len(intersection), len(B_only)), alpha=0.5, subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/total):1.0%}" + ")")
    #venn2_circles(subsets=(len(A_only), len(intersection), len(B_only)), linewidth=1.0)
    print(f'Saving to: {args.O}/{sample_A_name}_and_{sample_B_name}_percent.svg')
    # Save venn diagram to file
    plt.savefig(f'{args.O}/{sample_A_name}_and_{sample_B_name}_percent.svg')
    plt.close()
    
    return None 

if __name__ == '__main__':
    main()