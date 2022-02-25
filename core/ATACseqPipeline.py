##############
# ATAC-Seq Pipline for automated analysis 
# When specifying paths, never add the trailing '/'
#
#
#
##############

import pandas as pd 
import numpy as np 
import glob, os 

class Pipeline(): 
    def __init__(self, data_path, dry_run=False): 
        # Maybe I should check Python version and installations here!?
        self.data_path = data_path
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
        self.submission_path = data_path + '/submissions'
        if not os.path.exists(self.submission_path):
            os.makedirs(self.submission_path)
        self.dry_run = dry_run
        return None 
    
    def from_ssheet(self, ssheet_path): 
        '''
        Configure Pipeline from tab-delimited spreadsheet. Must contain the 
        following columns: SampleName, Read1 (Path to read 1 of FASTQ file), 
        Read2 (Path to read 2 of FASTQ file), Status (Atleast one sample
        must be "Control", TechRep (Technical Replicate number)
        '''
        self.ssheet = pd.read_csv(ssheet_path, sep='\t')
        print(self.ssheet.groupby('Status').count()['TechRep'])
        return None 

    def align_fastqs(self, genome_path):
        '''
        Align fastq files using Bowtie2. Genome must be pre-indexed. 
        #SBATCH -J bowtie2_job
        #SBATCH --output=slurm-outputs/slurm-%j-%x.out
        ''' 
        if not os.path.exists(self.data_path + '/bams/'):
            os.makedirs(self.data_path + '/bams/')
            print(self.data_path + '/bams/')
        
        for index, sample in self.ssheet.iterrows():
            sample_name = sample['SampleName'] + '.bam'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
module load bowtie2;
module load samtools;
bowtie2 -x {genome_path} -1 {sample['Read1']} -2 {sample['Read2']} | samtools sort -o {self.data_path + '/bams/' +  sample_name};
            '''
            f = open(f'{self.submission_path}/align_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/align_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/align_{sample["SampleName"]}.sh')
        return None

    def remove_duplicates(self):
        '''
        Remove duplicates that are an artefact of PCR. This must be done in a paired-end sequencing manner. 
        Picard's MarkDups is good for paired end sequencing reads. Much more "intelligent" than samtools rmdups. 
        '''
        if not os.path.exists(self.data_path + '/bams_noDups/'):
            os.makedirs(self.data_path + '/bams_noDups/')
            print(self.data_path + '/bams_noDups/')

        for bam in os.listdir(f'{self.data_path + "/bams/"}'):
            if bam.endswith('.bam'):
                cleaned_bam = f'{self.data_path + "/bams_noDups/" + bam}'
                metrics_file = f'{self.data_path + "/MarkDups_Metrics/" + bam + ".metrics"}'
                cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
module load gatk;
gatk MarkDuplicates -I {bam} -O {cleaned_bam} -M {metrics_file} --REMOVE_DUPLICATES=TRUE;        
                '''
        return None 
    
    def remove_mito(self): 
        '''
        Remove sequences that align to mitochondrial DNA. Usually overepressented because mitochondrial DNA is devoid of histones, 
        so Tn5 transposase is very active here.  
        '''
        if not os.path.exists(self.data_path + '/bams_noMito/'):
            os.makedirs(self.data_path + '/bams_noDups_noMito/')
            print(self.data_path + '/bams_noDups_noMito/')
        return None 

    def quality_control(self):
        '''
        Print FASTQC reports. 
        '''
        return None 

    