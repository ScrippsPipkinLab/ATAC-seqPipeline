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
import re 
import matplotlib.pyplot as plt
import shutil
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
        if os.path.exists(self.submission_path):
            shutil.rmtree(self.submission_path)
        os.makedirs(self.submission_path)     

        self.dry_run = dry_run
        self.output_file_dir = data_path + '/slurm_outputs'
        if not os.path.exists(self.output_file_dir):
            os.makedirs(self.output_file_dir)
        
        # !!!!!!!!!  Load the conda environment here
        return None 
    
    def configure_envs(self):
        os.system('conda env create --file envs/ATACseq_env.yml')
        return None 
    
    def from_ssheet(self, ssheet_path): 
        '''
        Configure Pipeline from tab-delimited spreadsheet. Must contain the 
        following columns: SampleName, Read1 (Path to read 1 of FASTQ file), 
        Read2 (Path to read 2 of FASTQ file), Status (Atleast one sample
        must be "Control", TechRep (Technical Replicate number)
        '''
        self.ssheet_path = ssheet_path
        self.ssheet = pd.read_csv(ssheet_path, sep='\t')

        return None 

    def align_fastqs(self, genome_path):
        '''
        Align fastq files using Bowtie2. Genome must be pre-indexed. 
        ''' 
        if not os.path.exists(self.data_path + '/bams/'):
            os.makedirs(self.data_path + '/bams/')
            print(self.data_path + '/bams/')
        
        for index, sample in self.ssheet.iterrows():
            sample_name = sample['SampleName'] + '.bam'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=256000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=24:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load bowtie2;
module load samtools;
bowtie2 -x {genome_path} -1 {sample['Read1']} -2 {sample['Read2']} | samtools sort -o {self.data_path + '/bams/' +  sample_name};
            '''
            f = open(f'{self.submission_path}/align_fastqs_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/align_fastqs_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/align_fastqs_{sample["SampleName"]}.sh')
        return None

    def legacy_remove_duplicates(self):
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
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
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
            bam = f'{self.data_path + "/bams/" + sample["SampleName"] + ".bam"}'
            bam_noMito = f'{self.data_path + "/bams_noMito/" + sample["SampleName"] + ".bam"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
python {self.app_path + '/bam_whitelist.py'} {bam} {whitelist_file} {bam_noMito};
            '''
            f = open(f'{self.submission_path}/remove_mito_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/remove_mito_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/remove_mito_{sample["SampleName"]}.sh')
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
            bam_list.append(bam)
        bam_list = ' '.join(bam_list)
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load fastqc;
fastqc {bam_list} -o {self.data_path + '/fastqc/'};
        '''
        f = open(f'{self.submission_path}/quality_control_{sample["SampleName"]}.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/quality_control_{sample["SampleName"]}.sh')
        else:
            os.system(f'sbatch {self.submission_path}/quality_control_{sample["SampleName"]}.sh')
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
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load samtools; 
samtools mpileup {bam} -o {pileup};            
            '''
            f = open(f'{self.submission_path}/pileup_reads_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/pileup_reads_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/pileup_reads_{sample["SampleName"]}.sh')
        return None 
    
    def call_peaks_joint(self): 
        '''
        Use MACS3 to call peaks. MACS3 now does a lot of the downstream processing (peak calling,
        differential peaks & basic plotting). This function uses MACS3 to call peaks on treatment 
        over control. Unsure if this is the best way to do this... keeping this until comparisons can be 
        made.
        '''
        if not os.path.exists(self.data_path + '/macs3/'):
            os.makedirs(self.data_path + '/macs3/')
        # Assign controls 
        c = np.unique((self.ssheet[self.ssheet['C/T']=='C']['Status']).to_numpy())
        t = np.unique((self.ssheet[self.ssheet['C/T']=='T']['Status']).to_numpy())
        self.experimental_combs = list(itertools.product(c, t))
        for control, treatment in self.experimental_combs: 
            control_files = (self.ssheet[self.ssheet['Status']==control]['SampleName']).to_numpy()
            control_files = [f'{self.data_path + "/bams_noMito/" + bam_file + ".bam"}' for bam_file in control_files]
            control_files = ' '.join(control_files)
            treatment_files = (self.ssheet[self.ssheet['Status']==treatment]['SampleName']).to_numpy()
            treatment_files = [f'{self.data_path + "/bams_noMito/" + bam_file + ".bam"}' for bam_file in treatment_files]
            treatment_files = ' '.join(treatment_files)
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
macs3 callpeak -t {treatment_files} \
-c {control_files} --outdir {self.data_path + '/macs3/'} \
--name {control}_{treatment} --format BAMPE -g mm --nomodel \
-q 0.01 --call-summits -B;
            '''
            f = open(f'{self.submission_path}/call_peaks_joint_{control}_{treatment}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/call_peaks_joint_{control}_{treatment}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/call_peaks_joint_{control}_{treatment}.sh')
        return None

    def remove_duplicates(self):
        '''
        Upgraded duplicate removal using samtools markdups. This plays better with fastq's that have optical duplicates. 
        '''
        list_of_dirs = ['bams_collated', 'bams_fixmated', 'bams_sorted', 'bams_noDups']
        for dir in list_of_dirs: 
            if not os.path.exists(self.data_path + f'/{dir}/'):
                os.makedirs(self.data_path + f'/{dir}/')
        
        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noMito/" + sample["SampleName"] + ".bam"}'
            collated = f'{self.data_path + "/bams_collated/" + sample["SampleName"] + ".bam"}'
            fixmated = f'{self.data_path + "/bams_fixmated/" + sample["SampleName"] + ".bam"}'
            sorted_ = f'{self.data_path + "/bams_sorted/" + sample["SampleName"] + ".bam"}'
            noDups = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'

            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=256000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load samtools; 
samtools collate --threads 16 -o {collated} {bam} int_dir{sample["SampleName"]};
samtools fixmate -m {collated} {fixmated};
samtools sort -o {sorted_} {fixmated};
samtools markdup -r {sorted_} {noDups};            
            '''
            f = open(f'{self.submission_path}/remove_duplicates_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/remove_duplicates_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/remove_duplicates_{sample["SampleName"]}.sh')
        return None

    def get_stats(self):
        '''
        Calculate statistics (% mapped, % duplicate, % Mitochondrial) from BAM
        '''
        if not os.path.exists(self.data_path + '/stats/'):
            os.makedirs(self.data_path + '/stats/')
        if not os.path.exists(self.data_path + '/flagstats/'):
            os.makedirs(self.data_path + '/flagstats/')
        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams/" + sample["SampleName"] + ".bam"}'
            noDups = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            noMito = f'{self.data_path + "/bams_noMito/" + sample["SampleName"] + ".bam"}'
            peaks = f'{self.data_path + "/macs3/" + sample["SampleName"] + "_peaks.narrowPeak"}'

            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load samtools;
# Total and mapped reads
samtools flagstat {bam} > {self.data_path + "/flagstats/" + sample["SampleName"] + ".flagstat"};
# Duplicate Reads 
samtools view -c {noDups} > {self.data_path + "/stats/" + sample["SampleName"] + ".dups"};
# Mitochondrial Reads
samtools view -c {noMito} > {self.data_path + "/stats/" + sample["SampleName"] + ".mito"};
# Number of Peaks 
wc -l {peaks} > {self.data_path + "/stats/" + sample["SampleName"] + ".peaks"};
            '''
            f = open(f'{self.submission_path}/get_stats_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/get_stats_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/get_stats_{sample["SampleName"]}.sh')
        return None 
    
    def legacy_call_peaks(self):
        '''
        Use MACS3 to call peaks. MACS3 now does a lot of the downstream processing (peak calling,
        differential peaks & basic plotting)
        '''
        if not os.path.exists(self.data_path + '/macs3/'):
            os.makedirs(self.data_path + '/macs3/')
        ##### NEW 
        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
macs3 callpeak -t {bam} \
--outdir {self.data_path + '/macs3/'} \
--name {sample["SampleName"]} --format BAMPE -g mm --nomodel \
-p 0.1 --call-summits -B --shift -100 --extsize 200;
            '''
            f = open(f'{self.submission_path}/call_peaks_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/call_peaks_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/call_peaks_{sample["SampleName"]}.sh') 
        return None 

    def call_peaks_macs2(self):
        '''
        Use MACS2 to call peaks. MACS2 now does a lot of the downstream processing (peak calling,
        differential peaks & basic plotting)
        Trying out Harvard guidelines
        '''
        if not os.path.exists(self.data_path + '/macs2/'):
            os.makedirs(self.data_path + '/macs2/')
        ##### NEW 
        for index, sample in self.ssheet.iterrows():
            bed = f'{self.data_path + "/shifted_peaks/" + sample["SampleName"] + ".bed"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load macs;
macs2 callpeak -t {bed} \
--outdir {self.data_path + '/macs2/'} \
--name {sample["SampleName"]} --format BED -g mm -B --shift -100 --extsize 200 --SPMR -p 0.01;
            '''
            f = open(f'{self.submission_path}/call_peaks_macs2_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/call_peaks_macs2_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/call_peaks_macs2_{sample["SampleName"]}.sh') 
        return None 

    def bam_to_bed(self):
        '''
        Convert BAM to BED
        '''
        if not os.path.exists(self.data_path + '/beds/'):
            os.makedirs(self.data_path + '/beds/')
        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            bed = f'{self.data_path + "/beds/" + sample["SampleName"] + ".bed"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load bedtools;
bedtools bamtobed -i {bam} > {bed};
            '''
            f = open(f'{self.submission_path}/bam_to_bed_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/bam_to_bed_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/bam_to_bed_{sample["SampleName"]}.sh')
        return None
            
    def fastqc(self):
        '''
        Run FastQC on BAM files. 
        '''
        if not os.path.exists(self.data_path + '/fastqc/'):
            os.makedirs(self.data_path + '/fastqc/')
        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
fastqc {bam} -o {self.data_path + '/fastqc/'};
            '''
            f = open(f'{self.submission_path}/fastqc_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/fastqc_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/fastqc_{sample["SampleName"]}.sh')

        return None 

    def multiqc_report(self):
        '''
        Generate a multiqc report
        ''' 
        if not os.path.exists(self.data_path + '/multiqc/'):
            os.makedirs(self.data_path + '/multiqc/')

        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
multiqc {self.data_path + '/fastqc/'} -o {self.data_path + '/multiqc/'};
            '''
        f = open(f'{self.submission_path}/multiqc.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/multiqc.sh')
        else:
            os.system(f'sbatch {self.submission_path}/multiqc.sh')
        
        return None

    def bed_to_gtf(self): 
        '''
        Convert merged BED peaks  to GTF for use in feature Counts 
        '''
        bed = f'{self.data_path + "/merged/merged_peaks_HOMER.bed"}'
        gtf = f'{self.data_path + "/merged/merged_peaks_HOMER.gtf"}'
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --output={self.output_file_dir}/%j.out
#SBATCH --time=1:00:00
source ~/.bashrc;
source activate ATACseq_env;
module load python3;
python3 {self.app_path}/bed2gtf.py -i {bed} -o {gtf};
        '''
        f = open(f'{self.submission_path}/bed_to_gtf.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/bed_to_gtf.sh')
        else:
            os.system(f'sbatch {self.submission_path}/bed_to_gtf.sh')
        return None 

    def peak_counts(self):
        '''
        Generate counts from bam using subread featureCounts.
        WAY faster than HTseq-count 
        '''
        if not os.path.exists(self.data_path + '/counts/'):
            os.makedirs(self.data_path + '/counts/')
        merged_peaks = f'{self.data_path + "/merged/merged_peaks_HOMER.gtf"}'
        bam_files = []
        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            bam_files.append(bam)
        bam_files = ' '.join(bam_files)
        counts = f'{self.data_path + "/counts/peak_counts.counts"}'
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load subread;
featureCounts -T 16 -a {merged_peaks} -t 'peak' -g 'peak_id' -p -o {counts} {bam_files};
        '''
        f = open(f'{self.submission_path}/peak_counts.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/peak_counts.sh')
        else:
            os.system(f'sbatch {self.submission_path}/peak_counts.sh')
        return None 
    
    
    def legacy_merge_peaks(self): 
        '''
        Merge Peaks to create a single peak file. "union" of peaks, not intersection. 
        '''
        if not os.path.exists(self.data_path + '/merged/'):
            os.makedirs(self.data_path + '/merged/')
        all_peaks = []
        for index, sample in self.ssheet.iterrows():
            peaks = f'{self.data_path + "/macs2/" + sample["SampleName"] + "_peaks.narrowPeak"}'
            all_peaks.append(peaks)
        all_peaks = ' '.join(all_peaks)
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load bedtools; 
cat {all_peaks} | sort -k1,1 -k2,2n | bedtools merge -i - > {self.data_path + "/merged/" + "merged_peaks.bed"};
        '''
        f = open(f'{self.submission_path}/merge_peaks.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/merge_peaks.sh')
        else:
            os.system(f'sbatch {self.submission_path}/merge_peaks.sh')
        return None 

    def merge_peaks(self): 
        '''
        Merge Peaks to create a single peak file. "union" of peaks, not intersection. use HOMER mergePeaks
        '''
        if not os.path.exists(self.data_path + '/merged/'):
            os.makedirs(self.data_path + '/merged/')
        all_peaks = []
        for index, sample in self.ssheet.iterrows():
            peaks = f'{self.data_path + "/no_blacklist/" + sample["SampleName"] + "_peaks.narrowPeak"}'
            all_peaks.append(peaks)
        all_peaks = ' '.join(all_peaks)
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load homer; 
mergePeaks -d 200 {all_peaks} > {self.data_path + "/merged/" + "merged_peaks_HOMER.bed"};
        '''
        f = open(f'{self.submission_path}/merge_peaks.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/merge_peaks.sh')
        else:
            os.system(f'sbatch {self.submission_path}/merge_peaks.sh')
        return None

    def filter_blacklist(self, blacklist_file): 
        '''
        Remove regions that are known to interfere with ATAC-seq & Chip-seq 
        experiments. 
        bedtools intersect -v -a -b 
        '''
        if not os.path.exists(self.data_path + '/no_blacklist/'):
            os.makedirs(self.data_path + '/no_blacklist/')
        for index, sample in self.ssheet.iterrows(): 
            peaks = f'{self.data_path + "/macs2/" + sample["SampleName"] + "_peaks.narrowPeak"}'
            filtered = f'{self.data_path + "/no_blacklist/" + sample["SampleName"] + "_peaks.narrowPeak"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load bedtools; 
bedtools intersect -v -a {peaks} -b {blacklist_file} > {filtered};
            '''
            f = open(f'{self.submission_path}/filter_blacklist_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/filter_blacklist_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/filter_blacklist_{sample["SampleName"]}.sh') 
        return None 

    def get_insert_sizes(self): 
        if not os.path.exists(self.data_path + '/insert_sizes/'):
            os.makedirs(self.data_path + '/insert_sizes/')
        for index, sample in self.ssheet.iterrows(): 
            bam = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            insert_sizes = f'{self.data_path + "/insert_sizes/" + sample["SampleName"] + ".sizes"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load samtools;
samtools view {bam} | cut -f9 > {insert_sizes};
            '''
            f = open(f'{self.submission_path}/get_insert_sizes_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/get_insert_sizes_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/get_insert_sizes_{sample["SampleName"]}.sh') 
        return None 
    
    def plot_insert_sizes(self):
        '''
        Read in insert sizes (list of integers) and plot density function.
        '''
        if not os.path.exists(self.data_path + '/insert_sizes/'):
            os.makedirs(self.data_path + '/insert_sizes/')
        if not os.path.exists(self.data_path + '/figs/'):
            os.makedirs(self.data_path + '/figs/')
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
python plot_insertsizes.py {self.data_path} {self.ssheet_path};
        '''
        f = open(f'{self.submission_path}/plot_insert_sizes.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/plot_insert_sizes.sh')
        else:
            os.system(f'sbatch {self.submission_path}/plot_insert_sizes.sh')
        return None 

#     def optimize_p_value_macs(self, start=0.01, end=0.25, num=25):
#         '''
#         *** RUN WITH CAUTION: HUGE RUNTIME INCREASE ***
#         Optimize p-value cutoff for MACS peak calling. Will search all p-values within a given
#         (start, stop, # of vals) range. Uses numpy.linspace to generate a list of q-values.
#         '''
#         if not os.path.exists(self.data_path + '/optimize_pval/'):
#             os.makedirs(self.data_path + '/optimize_pval/')
#         else:
#             shutil.rmtree(self.data_path + '/optimize_pval/')
#             os.makedirs(self.data_path + '/optimize_pval/')

#         self.pval_range = np.linspace(start, end, num)
#         for index, sample in self.ssheet.iterrows():
#             bed = f'{self.data_path + "/beds/" + sample["SampleName"] + ".bed"}'
#             for pval in self.pval_range:
#                 cmd = f'''#! /usr/bin/bash
# #SBATCH --cpus-per-task=8
# #SBATCH --mem=64000
# #SBATCH --output={self.data_path + '/optimize_pval/'}/%j_{sample["SampleName"]}_{pval}.out
# source ~/.bashrc;
# source activate ATACseq_env;
# macs2 callpeak -t {bed} \
# --outdir {self.data_path + '/optimize_pval/'} \
# --format BED -g mm -B --shift -100 -p {pval} \
# --extsize 200 --SPMR -n {sample["SampleName"]}_{pval};
#                 '''
#                 f = open(f'{self.submission_path}/optimize_p_value_macs_{sample["SampleName"] + "_" + "pval_" + str(pval)}.sh', 'w+')
#                 f.write(cmd)
#                 f.close()
#                 if self.dry_run:
#                     print(f'created {self.submission_path}/optimize_p_value_macs_{sample["SampleName"] + "_" + "pval_" + str(pval)}.sh')
#                 else:
#                     os.system(f'sbatch {self.submission_path}/optimize_p_value_macs_{sample["SampleName"] + "_" + "pval_" + str(pval)}.sh')
#         return None

    def shift_peaks(self, genome_sizes): 
        '''
        Shift peaks by +4 and -5 bp for +ve and -ve strand. This has to be done to account for the offset 
        of Tn5 transposase. 
        '''
        if not os.path.exists(self.data_path + '/shifted_peaks/'):
            os.makedirs(self.data_path + '/shifted_peaks/')
        for index, sample in self.ssheet.iterrows():
            bed = f'{self.data_path + "/beds/" + sample["SampleName"] + ".bed"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load bedtools; 
bedtools shift -i {bed} -g {genome_sizes} -p 4 -m -5 > {self.data_path + "/shifted_peaks/" + sample["SampleName"] + ".bed"};
            '''
            f = open(f'{self.submission_path}/shift_peaks_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/shift_peaks_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/shift_peaks_{sample["SampleName"]}.sh')
        return None 
    
    # def return_best_pval(self):
    #     '''
    #     Evaluate the best p-value and pick out the resultant files. Default to a p-value of 0.05
    #     if "elbow" cannot be established. 
    #     '''
    #     if not os.path.exists(self.data_path + '/optimize_pval/'):
    #         raise FileNotFoundError('Optimize_pval directory not found. Run optimize_p_value_macs() first.')

    #     if not os.path.exists(self.data_path + '/slurm_outputs/'):
    #         raise FileNotFoundError('Cannot find slurm outputs directory.')
        
    #     # loop through files and find best p-value
    #     for index, sample in self.ssheet.iterrows():
    #         fname = str(sample['SampleName'])
    #         for peakcall in glob.glob(f'{self.data_path + "/optimize_pval/*" + fname + "*.xls"}'):
    #             try: 
    #                 f = open(peakcall, 'r')
    #                 pval = re.match(r'# name = (\w+)_pval_(\d.+)').group(2)
    #                 mean_frag_size = re.match(r'mean fragment size is determined as (\d.+) bp from treatment').group(1)
    #                 f.close()
    #             except:
    #                 pass
    #     return None

    def bam_to_bigWig(self):
        '''
        Make bigWig files from bed files. BigWigs are used for visualization in IGV.
        '''
        if not os.path.exists(self.data_path + '/bigWigs/'):
            os.makedirs(self.data_path + '/bigWigs/')

        for index, sample in self.ssheet.iterrows():
            bam = f'{self.data_path + "/bams_noDups/" + sample["SampleName"] + ".bam"}'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load samtools;
module load deeptools;
samtools index {bam};
bamCoverage -b {bam} -o {self.data_path + "/bigWigs/" + sample["SampleName"] + ".bw"} \
--binSize 10 -p 8 --normalizeUsing RPKM -bs 10; 
            '''
            f = open(f'{self.submission_path}/bam_to_bigWig_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/bam_to_bigWig_{sample["SampleName"]}.sh')
            else:
                os.system(f'sbatch {self.submission_path}/bam_to_bigWig_{sample["SampleName"]}.sh')
        return None

    def differential_analysis(self, differential_exp_script_path):
        '''
        Find differentially accessible peaks using DESeq2
        '''
        if not os.path.exists(self.data_path + '/deseq2/'):
            os.makedirs(self.data_path + '/deseq2/')
        if differential_exp_script_path is None:
            differential_exp_script_path = self.app_path + '/differential_exp.r'
        
        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
#SBATCH --time=2:00:00
source ~/.bashrc;
source activate ATACseq_env;
module load R;
Rscript {differential_exp_script_path} {self.data_path + '/counts/peak_counts.counts'} {self.ssheet_path} {self.data_path + '/deseq2/'};
        '''
        f = open(f'{self.submission_path}/differential_analysis.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/differential_analysis.sh')
        else:
            os.system(f'sbatch {self.submission_path}/differential_analysis.sh')
        return None 

    def assign_peaks_to_genes(self, gene_annotation_path):
        '''
        cut -f 2,3,4,5,6 merged_peaks_HOMER.bed | tail -n +2 | sort -k1,1V -k2,2n -k3,3n > merged_peaks_HOMER.formatted.bed;
        bedtools closest -a merged_peaks_HOMER.formatted.bed -b /blue/m.pipkin/s.nagaraja/MusRef/Mouse39Genes.bed -D -t first > matched_genes.txt;
        '''
        if not os.path.exists(self.data_path + '/deseq2/'):
            try: 
                raise FileNotFoundError
            except Exception as error: 
                (repr(error))

        cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=10:00:00
#SBATCH --output={self.output_file_dir}/%j.out
#SBATCH --time=2:00:00
module load bedtools; 
cut -f 2,3,4,5,6 {self.data_path + '/merged/'}merged_peaks_HOMER.bed | tail -n +2 | sort -k1,1V -k2,2n -k3,3n > {self.data_path + '/merged/'}merged_peaks_HOMER.formatted.bed;
bedtools closest -a {self.data_path + '/merged/'}merged_peaks_HOMER.formatted.bed -b {gene_annotation_path} -D -t first > {self.data_path + '/merged/'}matched_genes.txt;
        '''
        f = open(f'{self.submission_path}/assign_peaks_to_genes.sh', 'w+')
        f.write(cmd)
        f.close()
        if self.dry_run:
            print(f'created {self.submission_path}/assign_peaks_to_genes.sh')
        else:
            os.system(f'sbatch {self.submission_path}/assign_peaks_to_genes.sh')

        return None 

    def draw_venn_diagrams(self):
        '''
        Draw venn diagrams from narrowPeak files. Draw venn for every permutation of narrowpeak files.
        '''
        if not os.path.exists(self.data_path + '/venn/'):
            os.makedirs(self.data_path + '/venn/')
            
        peaks = [f'{self.data_path + "/macs2/" + sample["SampleName"] + "_peaks.narrowPeak"}' for index, sample in self.ssheet.iterrows()]
        print(peaks)
        permutes = list(itertools.permutations(peaks, 2))
        print(len(permutes))

        for permute in permutes:
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=16000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --output={self.output_file_dir}/%j.out
#SBATCH --time=2:00:00
source activate ATACseq_env;
python {self.app_path}/venn_diagrams.py -A {permute[0]} -B {permute[1]} -O {self.data_path + '/venn'};
    '''
            f = open(f'{self.submission_path}/draw_venn_diagrams.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/draw_venn_diagrams.sh')
            else:
                os.system(f'sbatch {self.submission_path}/draw_venn_diagrams.sh')

        return None

    def main(self): 
        '''
        Create a dag of job dependencies and submit. 
        '''
        if not os.path.exists(self.data_path + '/figs/'):
            os.makedirs(self.data_path + '/figs/')

        #### FIXED IN MAIN
        pipe = ["align_fastqs",
                "remove_mito",
                "remove_duplicates",
                "bam_to_bed",
                "shift_peaks",
                "call_peaks_macs2",
                "filter_blacklist",
                "merge_peaks",
                "bed_to_gtf",
                "peak_counts",
                "bam_to_bigWig"]
        dependency = {}
        for idx, proc in enumerate(pipe):
            proc_jobs = []
            if idx == 0: 
                for submission in glob.glob(f'{self.submission_path}/{proc}*.sh'):
                    reply = os.popen(f'sbatch {submission}').read()
                    if re.match(r'Submitted batch job (\d+)', reply):
                            job_num = re.match(r'Submitted batch job (\d+)', reply).group(1)
                            proc_jobs.append(job_num)
                    else:
                        print(f'{submission} failed to submit')
            else: 
                for submission in glob.glob(f'{self.submission_path}/{proc}*.sh'):
                    dep_string = ':'.join(dependency[pipe[idx-1]])
                    reply = os.popen(f'sbatch --dependency=afterok:{dep_string} {submission}').read()
                    #print(reply)
                    if re.match(r'Submitted batch job (\d+)', reply):
                            job_num = re.match(r'Submitted batch job (\d+)', reply).group(1)
                            proc_jobs.append(job_num)
                    else:
                        print(f'{submission} failed to submit')
            dependency[proc] = proc_jobs
        print(dependency)        
        return None

class PipelineSE(Pipeline): 
    '''
    Alternate pipeline for single-ended data. 
    '''
    def align_fastqs(self, genome_path):
        '''
        Align fastq files using Bowtie2. Genome must be pre-indexed. 
        ''' 
        if not os.path.exists(self.data_path + '/bams/'):
            os.makedirs(self.data_path + '/bams/')
            print(self.data_path + '/bams/')
        
        for index, sample in self.ssheet.iterrows():
            sample_name = sample['SampleName'] + '.bam'
            cmd = f'''#! /usr/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=256000
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --time=24:00:00
#SBATCH --output={self.output_file_dir}/%j.out
source ~/.bashrc;
source activate ATACseq_env;
module load bowtie2;
module load samtools;
bowtie2 -x {genome_path} -U {sample['Read1']} | samtools sort -o {self.data_path + '/bams/' +  sample_name};
            '''
            f = open(f'{self.submission_path}/align_fastqs_{sample["SampleName"]}.sh', 'w+')
            f.write(cmd)
            f.close()
            if self.dry_run:
                print(f'created {self.submission_path}/align_fastqs_{sample["SampleName"]}.sh')
            else: 
                os.system(f'sbatch {self.submission_path}/align_fastqs_{sample["SampleName"]}.sh')
        return None