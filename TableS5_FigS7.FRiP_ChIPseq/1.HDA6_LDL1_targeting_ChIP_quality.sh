################################# Step 1. idr calculation
#/work1/home/joweihsieh/20190731_HDA6/ChIP-seq/20231207_data_quality

##### tapir idr.s
files1=(ChIP_hda6_1_peaks.narrowPeak ChIP_ldl1_1_peaks.narrowPeak)
files2=(ChIP_hda6_2_peaks.narrowPeak ChIP_ldl1_2_peaks.narrowPeak)

len=${#files1[@]}

process_files() {
    file1="$1"
    file2="$2"
    base1="${file1%_peaks.narrowPeak}"
    base2="${file2%_peaks.narrowPeak}"
    output_file="${base1}_${base2}_idr"
    output_log_file="${base1}_${base2}_idr.log"
    echo "Processing files: ${base1} and ${base2}"
    /usr/local/bin/idr --samples "$file1" "$file2" --input-file-type narrowPeak --rank p.value --output-file "$output_file" --plot --log-output-file "$output_log_file" --initial-sigma 2 --soft-idr-threshold 0.05 --peak-merge-method sum
    if [ $? -ne 0 ]; then
        echo "Error processing files: ${base1} and ${base2}"
    fi
}


for ((i = 0; i < len; i++))
do
    process_files "${files1[$i]}" "${files2[$i]}"
done
################################# Step 1.2. collect percentage of idr (done)
#!/bin/bash


for log_file in *_idr.log
do
    sample_name=$(basename "$log_file" | cut -d '_' -f 1-2)
    percentage=$(grep "Number of peaks passing IDR cutoff of 0.05" "$log_file" | awk '{print $NF}') 
    echo "Sample: $sample_name, Percentage: $percentage" >> idr_percent.txt
done



################################# Step 2 FRiP (in hal)

#!/bin/bash

files1=(*_peaks.narrowPeak)

for file1 in "${files1[@]}"
do
        base1="${file1%_peaks.narrowPeak}"
        
        output_file="${base1}.SAF"
        output_file2="${base1}-frip.txt"
        input_bam="${base1}.bam"
        output_log="${base1}_friplog"
        echo "Processing files: ${base1}"
        awk 'OFS="\t" {print $1"-"$2+1"-"$3, $1, $2+1, $3, "+"}' $file1 > $output_file 
        featureCounts -p -a $output_file -F SAF -o $output_file2 $input_bam >> $output_log 2>&1

        if [ $? -ne 0 ]; then
            echo "Error processing files: ${base1}"
        fi
done

################################# Step 6.2 collect percentage of FRiP

#!/bin/bash


for log_file in *_friplog
do
    sample_name=$(basename "$log_file" | cut -d '_' -f 1-2)
    percentage=$(grep -oE 'Successfully assigned reads : [0-9]+ \(([0-9]+\.[0-9]+)%\)' "$log_file" | awk '{print $NF}')
    echo "Sample: $sample_name, Percentage: $percentage" >> FRiP_percent.txt
done
