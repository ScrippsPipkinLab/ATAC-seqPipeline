#############
# Remove reads that are not in a provided "whitelist" of chromosomes. 
# Make sure your whitelist has the same echromosome nomenclature as the reference genome.
# I may want to save the statistics printed at the end here, so that I can plot them later.  
#############


import sys 
import pysam 

bam = sys.argv[1]
whitelist = open(sys.argv[2], 'r').read().split('\n')
output = sys.argv[3]

def main(bam_file, whitelist_path, outfile):
    samfile = pysam.AlignmentFile(bam, "rb")
    outfile = pysam.AlignmentFile(output, "wb", template=samfile)
    alignment_count = 0
    garbage_alignment_count = 0
    for read in samfile.fetch(until_eof=True):
        if read.reference_name in whitelist:
            outfile.write(read)
            alignment_count += 1
        else: 
            garbage_alignment_count +=1 
    print(f'Chromosomal Reads Fraction:     {garbage_alignment_count/(garbage_alignment_count+alignment_count)}')
    print(f'Total Reads:                    {garbage_alignment_count+alignment_count}')
    return None

if __name__ == "__main__":
    main(bam, whitelist, output)