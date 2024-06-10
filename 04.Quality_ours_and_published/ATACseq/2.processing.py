
##################
# trimming
##################

import pandas as pd
import numpy as np
import time
import os
import re
import shutil
#from tqdm.notebook import tqdm
from subprocess import Popen,PIPE,call


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
    #run_code = Popen(command_TrimGalore, shell=False, stdout=PIPE, stderr=PIPE, encoding="utf-8")
    run_code = Popen(command_TrimGalore, shell=False, stdout=PIPE, stderr=PIPE)

    stdout, stderr = run_code.communicate()
    
    print(fastq1 + "\t finish")
    # log("*Finish TrimGalore, trimmed.fq file has been created!")
    # Finish running TrimGalore trimmed.fq created.
    # Move the .fastq_trimming_report.txt to TrimReport folder
    # Delete .fastq file
    # Clean_FastQFile()
    # Move_TrimReport()


files = [
    "SRR6410823", "SRR6410824", "SRR5874657", "SRR5874658", "SRR12344672", "SRR12344673",
    "SRR8742423", "SRR8742424", "SRR12265344", "SRR12265345", "SRR4000468", "SRR4000469",
    "SRR7512044", "SRR7512045", "SRR5874660", "SRR5874661", "SRR5829242", "SRR5829243",
    "SRR4000476", "SRR4000477"
]

for file in files:
    print("RunTrimGalore('{}_1_de.fastq.gz', '{}_2_de.fastq.gz')".format(file, file))
    RunTrimGalore('{}_1_de.fastq.gz'.format(file), '{}_2_de.fastq.gz'.format(file))

