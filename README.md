# ATAC-seqPipeline


This ATACseq pipeline provides several tools to analyze ATAC-seq experiments. This pipeline is configured to automatically submit jobs using the SLURM environment. 

__Installation:__
* conda \
Install with `conda env create --file envs/ATACseq_env.yml`
* SLURM environment configured to submit jobs. 
* Singularity (optional) \
If you want, you can use a Singularity environment with the environment pre-built. Install with `singularity build envs/Singularity`
* Github

To install from github, you can clone the repository using: \
`git clone https://github.com/ScrippsPipkinLab/ATAC-seqPipeline.git`

__Usage:__

Sample metadata is entered into a simple tab-delimited text file. Make sure to include the below columns in the correct order. See core/exp122_ssheet.txt for a template. 

1. SampleName
Unique sample name. Must be different for each replicate of each sample. Set this to Sample_ReplicateNumber when in doubt. 
2. Read1 \
Path to R1 FASTQ file. Future versions will support single-end reads. 
3. Read2 \
Path to R2 FASTQ file.
4. Status \
This respresents the group of replicates that represent the sample. The status is the same accross replicates. For example, this could be a RNAmir construct, or an organ, or even a cell type. 
5. CT \
Control or treatment status denoted by a C or T respectively.

See core/exp122_ssheet.txt for a template. 
The provided example notebook file goes through how to access the functions in an interactive Jupyter session. In brief, you can import the python library and setup your experiment with: 
```
import ATACseqPipeline
myexp = ATACseqPipeline.Pipeline(data_path='path/to/empty/dir', dry_run=True, app_path='/ATACseqPipeline')
myexp.from_ssheet(ssheet_path='ATACseqPipeline/data/exp122_ssheet.txt')
```

Once the sample sheet has been configured, the entire pipline can me run with:

```
myexp.main()
```
To view your submitted jobs as they complete, run the following in the shell. You can also run this directly in the notebook using "!" before each command. 

```
squeue -u your_username
```
Abort the running jobs using `core/cancel_job.sh` \
This pipeline will submit several jobs that are dependent on each other. The runtime can reach several hours depending on the available nodes. You can also run individual parts of the pipeline, and the example python notebook provides an in-depth walkthrough of this. The output from every step is saved. This will cause quite a large amount of data to be saved (> 50 gigabytes for a 6 sample experiment). _It is up to the user to delete files saved in_ `/data`.

![alt text]([https://github.com/ScrippsPipkinLab/ATAC-seqPipeline/ATACseq_Pipeline.pdf](https://github.com/ScrippsPipkinLab/ATAC-seqPipeline/blob/858bd5e9e639f8f96c4bece51e00dcd3fd23764c/ATACseq_Pipeline.pdf))

Please feel free to raise issues on Github or shoot me an email: \
Shashank Nagaraja \
snagaraja@scripps.edu
