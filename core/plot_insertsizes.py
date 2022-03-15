##########
# Plot insert sizes of paired-end reads
##########
import sys 
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import random
import pandas as pd 
import numpy as np
import seaborn as sns

data_path = sys.argv[1]
ssheet_path = sys.argv[2]

def main(sample_sheet_path, data_path): 
    ssheet = pd.read_csv(ssheet_path, sep='\t')
    for index, sample in ssheet.iterrows(): 
        size = f'{data_path + "/insert_sizes/" + sample["SampleName"] + ".sizes"}'
        # Read in list of insert sizes and get kernel density estimate using sklearn
        with open(size, 'r') as f:
            insert_sizes = f.readlines()
            insert_sizes = np.array([abs(int(x)) for x in insert_sizes])
            n,x,_ = plt.hist(insert_sizes, bins=400, range=(1, 1000), density=True, label=sample['SampleName'])
            bin_centers = 0.5*(x[1:]+x[:-1])
            fig, ax = plt.subplots()
            ax.plot(bin_centers, n)
            ax.semilogy(base=2)
            ax.set_title(f'Insert Size Distribution: {sample["SampleName"]}')
            ax.set_xlabel('Insert Size (bp)')
            ax.set_ylabel('Density')
            plt.savefig(f'{data_path + "/figs/" + sample["SampleName"] + "_InsertSizes.png"}')
            plt.close()
        f.close()
    return None 

if __name__ == "__main__":
    main(ssheet_path, data_path)