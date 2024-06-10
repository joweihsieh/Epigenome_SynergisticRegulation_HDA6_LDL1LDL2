####rmDup
library(magrittr)
Filelist=list.files(pattern=".gz")
datalist=list()

for (ff in seq(1,length(Filelist),2)){
        File1=Filelist[ff]
        File2=Filelist[(ff+1)]
        File1_f=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('_rmDup.fastq.gz'))
        File2_f=strsplit(as.character(File2),split='.',2)[[1]][1] %>% paste0(c('_rmDup.fastq.gz'))
        File_log=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('_rmDup_log.txt'))
        code=paste("rmDupPE.pl", File1, File2, File1_f, File2_f, ">", File_log)
        datalist[[ff]]=code
}

dataset=do.call(rbind,datalist)
write.table(dataset,"WGBS_rmDupPE.sh",quote=F,sep=" ",row.names=F,col.names=F)





#### trimming
library(magrittr)
Filelist=list.files(pattern="_rmDup.fastq.gz")
datalist=list()

for (ff in seq(1,length(Filelist),2)){
        File1=Filelist[ff]
        File2=Filelist[(ff+1)]
        File1_p=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('.paired.fastq.gz'))
        File2_p=strsplit(as.character(File2),split='.',2)[[1]][1] %>% paste0(c('.paired.fastq.gz'))
        File1_un=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('.unpaired.fastq.gz'))
        File2_un=strsplit(as.character(File2),split='.',2)[[1]][1] %>% paste0(c('.unpaired.fastq.gz'))
        code=paste("java -jar trimmomatic-0.33.jar PE -phred33", File1, File2, File1_p, File1_un, File2_p, File2_un ,"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
        datalist[[ff]]=code
}

dataset=do.call(rbind,datalist)
write.table(dataset,"WGBS_trimming.sh",quote=F,sep=" ",row.names=F,col.names=F)

#### fastqc
library(magrittr)
Filelist=list.files(pattern="_rmDup.paired.fastq.gz")
datalist=list()
for (ff in 1:length(Filelist)){
        File1=Filelist[ff]
        code=paste("fastqc", File1)
        datalist[[ff]]=code
}

dataset=do.call(rbind,datalist)
write.table(dataset,"WGBS_fastqc.sh",quote=F,sep=" ",row.names=F,col.names=F)


#### R1 mapping


library(magrittr)
Filelist=list.files(pattern="_rmDup.paired.fastq.gz")
datalist=list()
for (ff in seq(1,length(Filelist),2)){
        File1=Filelist[ff]
        output1=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('.paired.bam'))
        output2=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('.paired.lambda.bam'))
        code=paste("/usr/local/bin/BSseeker2-master/bs_seeker2-align.py -i", File1, "-g /work1/home/joweihsieh/genome/Arabidopsis/TAIR10/TAIR10.fa --aligner=bowtie2 --temp_dir=./ --bt2-p 8 -d /work1/home/joweihsieh/genome/Arabidopsis/TAIR10/BS_index/ -o", output1)
        code2=paste("/usr/local/bin/BSseeker2-master/bs_seeker2-align.py -i", File1, "-g ~/genome/lambda/lambda.fa --aligner=bowtie2 --temp_dir=./ --bt2-p 8 -d ~/genome/lambda/BS_index -o", output2)
        datalist[[ff]]=rbind(code,code2)
}

dataset=do.call(rbind,datalist)
write.table(dataset,"WGBS_mapping.sh",quote=F,sep=" ",row.names=F,col.names=F)

#### R2 rev

#### rc R2 mapping

library(magrittr)
Filelist=list.files(pattern="_rmDup.paired.rev.fastq.gz")
datalist=list()
for (ff in seq(1:length(Filelist))){
        File1=Filelist[ff]
        output1=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('.paired.rev.bam'))
        output2=strsplit(as.character(File1),split='.',2)[[1]][1] %>% paste0(c('.paired.rev.lambda.bam'))
        code=paste("/usr/local/bin/BSseeker2-master/bs_seeker2-align.py -i", File1, "-g /work1/home/joweihsieh/genome/Arabidopsis/TAIR10/TAIR10.fa --aligner=bowtie2 --temp_dir=./ --bt2-p 8 -d /work1/home/joweihsieh/genome/Arabidopsis/TAIR10/BS_index/ -o", output1)
        code2=paste("/usr/local/bin/BSseeker2-master/bs_seeker2-align.py -i", File1, "-g ~/genome/lambda/lambda.fa --aligner=bowtie2 --temp_dir=./ --bt2-p 8 -d ~/genome/lambda/BS_index -o", output2)
        datalist[[ff]]=rbind(code,code2)
}

dataset=do.call(rbind,datalist)
write.table(dataset,"WGBS_mapping_rev_R2.sh",quote=F,sep=" ",row.names=F,col.names=F)


#### call methylation
library(magrittr)
Filelist=list.files(pattern=".bam")
datalist=list()
for (ff in seq(1:length(Filelist))){
        File1=Filelist[ff]
        code=paste("/usr/bin/python2 /work1/home/yenmr/tools/BSseeker2/bs_seeker2-call_methylation.py -i ",File1,"-x -d /work1/home/joweihsieh/genome/Arabidopsis/TAIR10/BS_index/TAIR10.fa_bowtie2")
        datalist[[ff]]=code
}

dataset=do.call(rbind,datalist)
write.table(dataset,"WGBS_callmeth.sh",quote=F,sep=" ",row.names=F,col.names=F)



