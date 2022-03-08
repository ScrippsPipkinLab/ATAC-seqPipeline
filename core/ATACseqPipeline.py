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
import itertools


class Pipeline(): 
    def __init__(self, data_path, app_path, dry_run=False, conda=''): 
        print(os.getcwd())
        # Maybe I should check Python version and installations here!?
        self.app_path = app_path
        self.data_path = data_path
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
        self.submission_path = data_path + '/submissions'
        if not os.path.exists(self.submission_path):
            os.makedirs(self.submission_path)
        self.dry_run = dry_run
        # !!!!!!!!!  Load the conda environment here
        return None 
    
    def from_ssheet(self, ssheet_path): 
        '''
        Configure Pipeline from tab-delimited spreadsheet. Must contain the 
        following columns: SampleName, Read1 (Path to read 1 of FASTQ file), 
        Read2 (Path to read 2 of FASTQ file), Status (Atleast one sample
        must be "Control", TechRep (Technical Replicate number)
        '''
        self.ssheet = pd.read_csv(ssheet_path, sep='\t')
        #####
        # print(self.ssheet.groupby('Status').count()['TechRep'])
        #####
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
        if not os.path.exists(self.data_path + '/MarkDups_Metrics/'):
            os.makedirs(self.data_path + '/MarkDups_Metrics/')
            print(self.data_path + '/MarkDups_Metrics/')

        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams/" + sample["SampleName"] + ".bam"}'
            bam_noDups = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            metrics_file = f'{self.data_path + "/MarkDups_Metrics/" + sample["SampleName"] + ".metrics"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
module load gatk;
gatk MarkDuplicates -I {bam} -O {bam_noDups} -M {metrics_file} --REMOVE_DUPLICATES TRUE;        
                '''
            f = open(f'{self.submission_path}/removeDups_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/removeDups_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/removeDups_{sample["SampleName"]}.sh')
        return None 
    
    def remove_mito(self, whitelist_file='RefSeq-whitelist.txt'): 
        '''
        Remove sequences that align to mitochondrial DNA. Usually overepressented because mitochondrial DNA is devoid of histones, 
        so Tn5 transposase is very active here.  
        '''

        if not os.path.exists(self.data_path + '/bams_noMito/'):
            os.makedirs(self.data_path + '/bams_noMito/')
            print(self.data_path + '/bams_noDups_noMito/')
        
        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            bam_noMito = f'{self.data_path + "/bams_noMito/" + sample["SampleName"] + ".bam"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
# conda activate ATACseq_env;
python {self.app_path + '/bam_whitelist.py'} {bam} {whitelist_file} {bam_noMito};
            '''
            f = open(f'{self.submission_path}/removeMito_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/removeMito_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/removeMito_{sample["SampleName"]}.sh')
        return None 

    def quality_control(self):
        '''
        Print FASTQC reports. 
        '''
        if not os.path.exists(self.data_path + '/fastqc/'):
            os.makedirs(self.data_path + '/fastqc/')

        bam_list = []
        for index, sample in self.ssheet.iterrows(): 
            bam = f'{self.data_path + "/bams_noMito/" + sample["SampleName"] + ".bam"}'
            bam_list = bam_list.append(bam)
        bam_list = ' '.join(bam_list)
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
module load fastqc;
fastqc {bam_list} -o {self.data_path + '/fastqc/'};
        '''
        f = open(f'{self.submission_path}/fastqc_{sample["SampleName"]}.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/fastqc_{sample["SampleName"]}.sh')
        else:
            os.system(f'sbatch {self.submission_path}/fastqc_{sample["SampleName"]}.sh')
        return None 
    
    def pileup_reads(self): 
        '''
        Pileup reads using samtools' mpileup. We could also consider implementing GATK's pileup in the future here. 
        '''
        if not os.path.exists(self.data_path + '/pileups/'):
            os.makedirs(self.data_path + '/pileups/')
            print(self.data_path + '/pileups/')

        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noMito/" + sample["SampleName"] + ".bam"}'
            pileup = f'{self.data_path + "/pileups/" + sample["SampleName"] + ".pileup"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
module load samtools; 
samtools mpileup {bam} -o {pileup};            
            '''
            f = open(f'{self.submission_path}/pileup_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/pileup_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/pileup_{sample["SampleName"]}.sh')
        return None 
    
    def call_peaks(self): 
        '''
        Use MACS3 to call peaks. MACS3 now does a lot of the downstream processing (peak calling,
        differential peaks & basic plotting)
        '''
        if not os.path.exists(self.data_path + '/macs3/'):
            os.makedirs(self.data_path + '/macs3/')
            print(self.data_path + '/macs3/')
        # Assign controls 
        c = np.unique((self.ssheet[self.ssheet['C/T']=='C']['Status']).to_numpy())
        t = np.unique((self.ssheet[self.ssheet['C/T']=='T']['Status']).to_numpy())
        combs = list(itertools.product(c, t))
        for control, treatment in combs: 
            control_files = (self.ssheet[self.ssheet['Status']==control]['SampleName']).to_numpy()
            control_files = [f'{self.data_path + "/bams_noMito/" + bam_file + ".bam"}' for bam_file in control_files]
            control_files = ' '.join(control_files)
            treatment_files = (self.ssheet[self.ssheet['Status']==treatment]['SampleName']).to_numpy()
            treatment_files = [f'{self.data_path + "/bams_noMito/" + bam_file + ".bam"}' for bam_file in treatment_files]
            treatment_files = ' '.join(treatment_files)
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
/gpfs/home/snagaraja/miniconda3/envs/ATACseq_env/bin/macs3 callpeak -t {treatment_files} -c {control_files};
            '''
            f = open(f'{self.submission_path}/pileup_{control}_{treatment}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/pileup_{control}_{treatment}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/pileup_{control}_{treatment}.sh')
        return None