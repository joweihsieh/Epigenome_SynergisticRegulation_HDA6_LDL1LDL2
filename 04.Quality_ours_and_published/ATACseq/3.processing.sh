# conda activate /work4/home/guanjun/anaconda3

##################
# index for mapping
##################
 
nohup bowtie2-build --threads 20 TAIR10.fa TAIR10 &

##################
# mapping
##################
 
mkdir -p 3.BAM

######
input_file="acc_list.txt"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

while IFS= read -r file1; do
    input_file1="${file1}_1_de_val_1.fq.gz"
    input_file2="${file1}_2_de_val_2.fq.gz"

    output_file1="${file1}.bam"
    output_file2="Step2.bowtie2mapping.log"
    echo "Mapping to the genome for ${file1}" >> "$output_file2"
    # 使用bowtie2比对
    if bowtie2 -p 20 -x /work4/home/joweihsieh/genome/Arabidopsis/TAIR10/TAIR10 -1 "./2.Trim_Fastq/$input_file1" -2 "./2.Trim_Fastq/$input_file2" 2>>"$output_file2" | samtools sort -o "./3.BAM/$output_file1"; then
        echo "File processing successful: $file1" >> "$output_file2"
    else
        echo "File processing failed: $file1" >> "$output_file2"
    fi
done < "$input_file"


##################
# indexing
##################

#files=(*.bam)   
for file in ./*.bam; do

    samtools index $file
done

#/work3/guanjun/JW
#samtools index -@ 20 SRR12265345.bam


##################
# remove Pt and Mt
##################

out_dir="./4.rmPtMt"

if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

input_file="acc_list.txt"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

while IFS= read -r file1; do
    input_file1="./3.BAM/${file1}.bam"
    output_file1="./4.rmPtMt/${file1}_rmMtPt.bam"
    output_file2="Step3.rmMtPt.log"
    
    echo "Mapping to the genome for ${file1}" >> "$output_file2"
    
    if samtools index "$input_file1"; then
        echo "Indexing successful for $file1" >> "$output_file2"
    else
        echo "Indexing failed for $file1" >> "$output_file2"
        exit 1
    fi
    
    if python /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py "$input_file1" "$output_file1" Mt,Pt 2>>"$output_file2"; then
        echo "File processing successful: $file1" >> "$output_file2"
    else
        echo "File processing failed: $file1" >> "$output_file2"
    fi
done < "$input_file"



#python /work4/home/guanjun/tools/ATACgraph2/ATACgraph2/python3/00_rmChr.py SRR12265344.bam SRR12265344.rmMtPt.bam Mt,Pt
#conda activate /work4/home/guanjun/anaconda3
#python3 /work4/home/joweihsieh/bin/ATACgraph2/python3/00_rmChr.py SRR12265345.bam SRR12265345_rmMt.bam Mt,Pt > SRR12265345_rmMt.log &
#19433

##################
# indexing
##################

#files=(*.bam)   
for file in ./*_rmM.bam; do

    samtools index -@ 20 $file
done


#samtools index -@ 20 SRR12265345_rmMt.bam

##################
# size selection
##################

out_dir="./5.Select_size"

if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

input_file="acc_list.txt"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

# 逐行处理输入文件中的文件名
while IFS= read -r file1; do
    input_file1="./4.rmPtMt/${file1}_rmM.bam"
    output_file1="./5.Select_size/${file1}_rmM_150bp.bam"
    output_file2="Step4.size.selection.log"
    
    echo "Size selection for ${file1}" >> "$output_file2"
    
    
    # 使用python脚本进行操作
    if python3 /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py "$input_file1" "$output_file1" 2>>"$output_file2"; then
        echo "File processing successful: $file1" >> "$output_file2"
    else
        echo "File processing failed: $file1" >> "$output_file2" 
    fi


done < "$input_file"

#python3 /work4/home/joweihsieh/bin/ATACgraph2/python3/02_selectFragSize.py ../4.rmPtMt/SRR12265345_rmMt.bam SRR12265345_rmM_150bp.bam 2>>../Step4.size.selection.log


##################
# annotation file
##################

cp /work1/home/yenmr/genome/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Genes/genes.gtf ./
mv genes.gtf TAIR10.gtf
python3 /work4/home/joweihsieh/bin/ATACgraph2/python3/02_gtftoBed.py TAIR10.gtf TAIR10


##################
# peak calling
##################

out_dir="./6.Peak"

if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

input_file="acc_list.txt"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

while IFS= read -r file1; do
    input_file1="./5.Select_size/${file1}_rmM_150bp.bam"
    annotation_file="TAIR10_gene_body_merge.bed"
    output="./6.Peak/${file1}"
    output_file2="Step5.callpeak.log"
    
    echo "peak calling for ${file1}" >> "$output_file2"
    
   
    if python3 /work4/home/joweihsieh/bin/ATACgraph2/python3/03_callPeak.py -p  /work1/home/hsieh/anaconda3/envs/peakcall/bin/ -s 1 "$input_file1" "$output" "$annotation_file" 2>>"$output_file2"; then
        echo "File processing successful: $file1" >> "$output_file2"
    else
        echo "File processing failed: $file1" >> "$output_file2"
    fi
done < "$input_file"


#python /work4/home/joweihsieh/bin/ATACgraph/script/03_callPeak.py -s 1 -c S27_sort_rmMtPt_L150bp.bam S23_sort_rmMtPt_L150bp.bam S23_peak ../annotation/TAIR10_gene_promoter_bed6.bed


##################
# convert to bw
##################


out_dir="./6.Peak"

if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

input_file="acc_list.txt"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

while IFS= read -r file1; do
    input_file1="./5.Select_size/${file1}_rmM_150bp.bam"
    output_file1="./5.Select_size/${file1}_rmM_150bp.bw"
    output="./6.Peak/${file1}"
    output_file2="Step5_1.bw.log"
    
    echo "bam to bw for ${file1}" >> "$output_file2"


    if bamCoverage -b "$input_file1" -o "$output_file1" --normalizeUsing RPKM --binSize 10; then
            echo "BAM coverage successful: $file1" >> "$output_file2"
    else
            echo "BAM coverage failed: $file1" >> "$output_file2"
    fi
done < "$input_file"


#python /work4/home/joweihsieh/bin/ATACgraph/script/03_genePlot.py S23_peak_peaks.narrowPeak S23_peak_coverage.bw ./annotation/TAIR10


##################
# QC - idr
##################

mkdir -p 8.QCindex
cd /work4/home/joweihsieh/20240424_At_meth_ATAC_pub/ATAC_raw/8.QCindex
cp /work4/home/joweihsieh/20211112_HDA6_ATAC/4.Select/*_sort_rmMtPt_150bp.bam ./
ln -s /work4/home/joweihsieh/20211112_HDA6_ATAC/4.Select/*_peak_peaks.narrowPeak ./


ln -s /work4/home/joweihsieh/20240424_At_meth_ATAC_pub/ATAC_raw/6.Peak/*peaks.narrowPeak ./
ln -s /work4/home/joweihsieh/20240424_At_meth_ATAC_pub/ATAC_raw/5.Select_size/*_rmM_150bp.bam ./

#conda activate ATAC
#conda install idr
#vi /work4/home/joweihsieh/miniconda3/envs/ATAC/lib/python3.8/site-packages/idr/idr.py
#Replace numpy.int with int
#/work4/home/joweihsieh/20240424_At_meth_ATAC_pub/ATAC_raw/8.QCindex
##### tapir idr.s

files1=(WT_2_peak_peaks.narrowPeak HDA6_1_peak_peaks.narrowPeak LDL_3_peak_peaks.narrowPeak HL_2_peak_peaks.narrowPeak 
    SRR12265344_peaks.narrowPeak SRR12344672_peaks.narrowPeak SRR4000468_peaks.narrowPeak SRR4000476_peaks.narrowPeak
    SRR5829242_peaks.narrowPeak SRR5874657_peaks.narrowPeak SRR5874660_peaks.narrowPeak SRR6410823_peaks.narrowPeak
    SRR7512044_peaks.narrowPeak SRR8742423_peaks.narrowPeak)

files2=(WT_3_peak_peaks.narrowPeak HDA6_4_peak_peaks.narrowPeak LDL_4_peak_peaks.narrowPeak HL_3_peak_peaks.narrowPeak
    SRR12265345_peaks.narrowPeak SRR12344673_peaks.narrowPeak SRR4000469_peaks.narrowPeak SRR4000477_peaks.narrowPeak
    SRR5829243_peaks.narrowPeak SRR5874658_peaks.narrowPeak SRR5874661_peaks.narrowPeak SRR6410824_peaks.narrowPeak
    SRR7512045_peaks.narrowPeak SRR8742424_peaks.narrowPeak)

len=${#files1[@]}

process_files() {
    file1="$1"
    file2="$2"
    base1="${file1%_peaks.narrowPeak}"
    base2="${file2%_peaks.narrowPeak}"
    output_file="${base1}_${base2}_idr"
    output_log_file="${base1}_${base2}_idr.log"
    echo "Processing files: ${base1} and ${base2}"
    idr --samples "$file1" "$file2" --input-file-type narrowPeak --rank p.value --output-file "$output_file" --plot --log-output-file "$output_log_file" --initial-sigma 2 --soft-idr-threshold 0.05 --peak-merge-method sum
    if [ $? -ne 0 ]; then
        echo "Error processing files: ${base1} and ${base2}"
    fi
}


for ((i = 0; i < len; i++))
do
    process_files "${files1[$i]}" "${files2[$i]}"
done


#idr --samples WT_2_peak_peaks.narrowPeak WT_3_peak_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file WT_2_WT_3_idr --plot --log-output-file WT_2_WT_3_idr.log --initial-sigma 2 --soft-idr-threshold 0.05 --peak-merge-method sum



#### collect percentage of idr (done)
#!/bin/bash


for log_file in *_idr.log
do
    sample_name=$(basename "$log_file" | cut -d '_' -f 1-2)
    percentage=$(grep "Number of peaks passing IDR cutoff of 0.05" "$log_file" | awk '{print $NF}') 
    echo "Sample: $sample_name, Percentage: $percentage" >> idr_percent.txt
done


##################
# QC - FRiP
##################

files1=(*_peaks.narrowPeak)

for file1 in "${files1[@]}"
do
        base1="${file1%_peaks.narrowPeak}"
        
        output_file="${base1}.SAF"
        output_file2="${base1}-frip.txt"
        input_bam="${base1}_rmM_150bp.bam"
        output_log="${base1}_friplog"
        echo "Processing files: ${base1}"
        awk 'OFS="\t" {print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}' $file1 > $output_file 
        featureCounts -p -a $output_file -F SAF -o $output_file2 $input_bam >> $output_log 2>&1

        if [ $? -ne 0 ]; then
            echo "Error processing files: ${base1}"
        fi
done

##### some files were previously named with sorted
#!/bin/bash

files1=(*peak_peaks.narrowPeak)

for file1 in "${files1[@]}"
do
        base1="${file1%_peak_peaks.narrowPeak}"
        
        output_file="${base1}.SAF"
        output_file2="${base1}-frip.txt"
        input_bam="${base1}_sort_rmMtPt_150bp.bam"
        output_log="${base1}_friplog"
        echo "Processing files: ${base1}"
        awk 'OFS="\t" {print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}' $file1 > $output_file 
        featureCounts -p -a $output_file -F SAF -o $output_file2 $input_bam >> $output_log 2>&1

        if [ $? -ne 0 ]; then
            echo "Error processing files: ${base1}"
        fi
done


#### collect percentage of FRiP

#!/bin/bash


for log_file in *_friplog
do
    sample_name=$(basename "$log_file" | cut -d '_' -f 1-2)
    percentage=$(grep 'Successfully assigned fragments :' "$log_file" | awk '{print $(NF-1)}')

    echo "Sample: $sample_name, Percentage: $percentage" >> FRiP_percent.txt
done



