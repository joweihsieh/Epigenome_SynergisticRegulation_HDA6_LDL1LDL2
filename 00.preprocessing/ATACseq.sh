# /Users/joweihsieh/Dropbox/Prof.Chen/project/ATAC/ATAC_RUN.ipynb
# pip3 install jupyter
# pip3 install ipython notebook
#jupyter notebook ATAC_RUN.ipynb

# md5sum *fq.gz  >> md5.txt


import pandas as pd
import numpy as np
import time
import os
import re
import shutil
from tqdm.notebook import tqdm
from subprocess import Popen,PIPE,call

###########################
# Tools' path
###########################
global INFO

INFO = {
    # location Infromation of the programs
    "parallel-fastq-dump":"parallel-fastq-dump",
    "trim-galore":"/work4/home/joweihsieh/bin/TrimGalore-0.6.10/trim_galore",
    "bowtie2":"bowtie2",
    "samtools":"/work1/home/hsieh/softwares/samtools-1.10/samtools",
    "macs2":"macs2",
    "bamCoverage":"bamCoverage",
    "GenomeIndex":"/work4/home/joweihsieh/genome/Arabidopsis/TAIR10/TAIR10"
}

###########################
#  Step 1
###########################

def RunTrimGalore(fastq1, fastq2):
    # log("*Start running TrimGalore program")
    
    ####################
    ###  Parameters  ###
    ####################
    
    program = INFO["trim-galore"]
    quality = 30
    cores = 30  # Exchange the threads for increasing the speed
    outdir = "./"
    fastq1 = "./1.Raw_Fastq/" + fastq1
    fastq2 = "./1.Raw_Fastq/" + fastq2
    # Run the codes
    
    command_TrimGalore = [program, "-q", str(quality), "--output_dir", outdir, "-j", str(cores), "--paired", fastq1,
                           fastq2, "--fastqc"]
    run_code = Popen(command_TrimGalore, shell=False, stdout=PIPE, stderr=PIPE)

    stdout, stderr = run_code.communicate()
    
    print(fastq1 + "\t finish")



files = [file for file in os.listdir("./1.Raw_Fastq") if re.search(r'.fq$', file)]
files


RunTrimGalore('ATAC_WT2_DKDL210010287-1a_HH7KFCCX2_L5_1.fq','ATAC_WT2_DKDL210010287-1a_HH7KFCCX2_L5_2.fq')
RunTrimGalore('ATAC_WT3_DKDL210010288-1a_HH7KFCCX2_L5_1.fq','ATAC_WT3_DKDL210010288-1a_HH7KFCCX2_L5_2.fq')


RunTrimGalore('ATAC_HDA6_1_DKDL210010289-1a_HH7KFCCX2_L5_1.fq','ATAC_HDA6_1_DKDL210010289-1a_HH7KFCCX2_L5_2.fq')
RunTrimGalore('ATAC_HDA6_4_DKDL210010290-1a_HH7KFCCX2_L5_1.fq','ATAC_HDA6_4_DKDL210010290-1a_HH7KFCCX2_L5_2.fq')


RunTrimGalore('ATAC_LDL_3_DKDL210010291-1a_HH7KFCCX2_L5_1.fq','ATAC_LDL_3_DKDL210010291-1a_HH7KFCCX2_L5_2.fq')
RunTrimGalore('ATAC_LDL_4_DKDL210010292-1a_HH7KFCCX2_L5_1.fq','ATAC_LDL_4_DKDL210010292-1a_HH7KFCCX2_L5_2.fq')


RunTrimGalore('ATAC_HL_2_DKDL210010293-1a_HH7KFCCX2_L5_1.fq','ATAC_HL_2_DKDL210010293-1a_HH7KFCCX2_L5_2.fq')
RunTrimGalore('ATAC_HL_3_DKDL210010294-1a_HH7KFCCX2_L5_1.fq','ATAC_HL_3_DKDL210010294-1a_HH7KFCCX2_L5_2.fq')


RunTrimGalore('gDNA_DKDL210010295-1a_HH7KFCCX2_L5_1.fq','gDNA_DKDL210010295-1a_HH7KFCCX2_L5_2.fq')

###########################
#  Step 2
###########################

#mkdir 2.Trim_Fastq/
#mv *val*.fq ./2.Trim_Fastq


def RunBowtie2(file_1,file_2,name):
#     trimfastq = GSM+"_trimmed.fq"
#     log("*Start read mapping by bowtie2 program")
    file_1 = "./2.Trim_Fastq/" + file_1
    file_2 = "./2.Trim_Fastq/" + file_2
    program = INFO["bowtie2"]
    genome = INFO["GenomeIndex"]
    threads = 30 # Exchange the threads for increasing the speed
        
    command_Bowtie2 = [program,"-x",genome,"-1",file_1,"-2",file_2,"-p",str(threads),"-S",name+".sam","2>"+name+".bowtie_log"]
    run_code = Popen(command_Bowtie2,shell=False,stdout=PIPE,stderr=PIPE,encoding="ascii")
    stdout, stderr = run_code.communicate()



files = [file for file in os.listdir("./2.Trim_Fastq") if re.search(r'.fq$', file)]
files



RunBowtie2("ATAC_WT2_DKDL210010287-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_WT2_DKDL210010287-1a_HH7KFCCX2_L5_2_val_2.fq","WT_2")
RunBowtie2("ATAC_WT3_DKDL210010288-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_WT3_DKDL210010288-1a_HH7KFCCX2_L5_2_val_2.fq","WT_3")


RunBowtie2("ATAC_HDA6_1_DKDL210010289-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_HDA6_1_DKDL210010289-1a_HH7KFCCX2_L5_2_val_2.fq","WT_2")
RunBowtie2("ATAC_HDA6_4_DKDL210010290-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_HDA6_4_DKDL210010290-1a_HH7KFCCX2_L5_2_val_2.fq","HDA6_4")


RunBowtie2("ATAC_LDL_3_DKDL210010291-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_LDL_3_DKDL210010291-1a_HH7KFCCX2_L5_2_val_2.fq","LDL_3")
RunBowtie2("ATAC_LDL_4_DKDL210010292-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_LDL_4_DKDL210010292-1a_HH7KFCCX2_L5_2_val_2.fq","LDL_4")


RunBowtie2("ATAC_HL_2_DKDL210010293-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_HL_2_DKDL210010293-1a_HH7KFCCX2_L5_2_val_2.fq","HL_2")
RunBowtie2("ATAC_HL_3_DKDL210010294-1a_HH7KFCCX2_L5_1_val_1.fq","ATAC_HL_3_DKDL210010294-1a_HH7KFCCX2_L5_2_val_2.fq","HL_3")

RunBowtie2("gDNA_DKDL210010295-1a_HH7KFCCX2_L5_1_val_1.fq","gDNA_DKDL210010295-1a_HH7KFCCX2_L5_2_val_2.fq","gDNA")


###########################
#  Step 3
###########################

def ConvertSAMtoBAM(file):

    program = INFO["samtools"]
    threads = 20 # Exchange the threads for increasing the speed
    os.system(program+" view -bS "+ file+".sam" +" -@ "+str(threads)+" > "+file+".bam")
    os.system("nohup "+ program +" sort -o "+file+"_sort.bam -@ "+str(threads)+" "+ file+".bam" +" >/dev/null 2>&1")
    os.system(program +" index "+ file +"_sort.bam")


ConvertSAMtoBAM("WT_2")
ConvertSAMtoBAM("WT_3")

ConvertSAMtoBAM("HDA6_1")
ConvertSAMtoBAM("HDA6_4")

ConvertSAMtoBAM("LDL_3")
ConvertSAMtoBAM("LDL_4")

ConvertSAMtoBAM("HL_2")
ConvertSAMtoBAM("HL_3")

ConvertSAMtoBAM("gDNA")

###########################
#  Step 4 - remove Pt and Mt
###########################

python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py WT_2_sort.bam WT_2_sort_rmMtPt.bam Mt,Pt
python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py WT_3_sort.bam WT_3_sort_rmMtPt.bam Mt,Pt

python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py HDA6_1_sort.bam HDA6_1_sort_rmMtPt.bam Mt,Pt
python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py HDA6_4_sort.bam HDA6_4_sort_rmMtPt.bam Mt,Pt

python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py LDL_3_sort.bam LDL_3_sort_rmMtPt.bam Mt,Pt
python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py LDL_4_sort.bam LDL_4_sort_rmMtPt.bam Mt,Pt

python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py HL_2_sort.bam HL_2_sort_rmMtPt.bam Mt,Pt
python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py HL_3_sort.bam HL_3_sort_rmMtPt.bam Mt,Pt

python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py gDNA_sort.bam gDNA_sort_rmMtPt.bam Mt,Pt

###########################
#  Step 5 - select 150bp
###########################

python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/raw_data/ATAC_WT2/WT_2_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/raw_data/ATAC_WT2/WT_2_sort_rmMtPt_150bp.bam
python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/new/raw_data/ATAC_WT3/WT_3_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/new/raw_data/ATAC_WT3/WT_3_sort_rmMtPt_150bp.bam

python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/raw_data/ATAC_LDL_3/LDL_3_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/raw_data/ATAC_LDL_3/LDL_3_sort_rmMtPt_150bp.bam
python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/raw_data/ATAC_LDL_4/LDL_4_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/raw_data/ATAC_LDL_4/LDL_4_sort_rmMtPt_150bp.bam


python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/raw_data/ATAC_HL_2/HL_2_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/raw_data/ATAC_HL_2/HL_2_sort_rmMtPt_150bp.bam
python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/raw_data/ATAC_HL_3/HL_3_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/raw_data/ATAC_HL_3/HL_3_sort_rmMtPt_150bp.bam


python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/raw_data/ATAC_HDA6_1/HDA6_1_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/raw_data/ATAC_HDA6_1/HDA6_1_sort_rmMtPt_150bp.bam
python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/new/raw_data/ATAC_HDA6_4/HDA6_4_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/new/raw_data/ATAC_HDA6_4/HDA6_4_sort_rmMtPt_150bp.bam


python /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ~/20211112_HDA6_ATAC/new/raw_data/gDNA/gDNA_sort_rmMtPt.bam ~/20211112_HDA6_ATAC/new/raw_data/gDNA/gDNA_sort_rmMtPt_150bp.bam


# indexing
files = [file for file in os.listdir("./") if re.search(r'_150bp.bam$', file)]
for file in files:
    os.system(INFO["samtools"] +" index "+ file)


###########################
#  Step 6 - peak calling
###########################


python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 WT_2_sort_rmMtPt_150bp.bam WT_2_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 WT_3_sort_rmMtPt_150bp.bam WT_3_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 LDL_3_sort_rmMtPt_150bp.bam LDL_3_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 LDL_4_sort_rmMtPt_150bp.bam LDL_4_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 HL_2_sort_rmMtPt_150bp.bam HL_2_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 HL_3_sort_rmMtPt_150bp.bam HL_3_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 HDA6_1_sort_rmMtPt_150bp.bam HDA6_1_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 HDA6_4_sort_rmMtPt_150bp.bam HDA6_4_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed
python /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -s 1 gDNA_sort_rmMtPt_150bp.bam gDNA_peak /work1/home/joweihsieh/20200813_ATAC-seq/first/annotation/TAIR10_gene_promoter_bed6.bed

