
##################
# Downloading
##################
# Read SRA accession numbers from acc_list.txt into an array

mapfile -t files1 < acc_list.txt
printf '%s\n' "${files1[@]}"

process_files() {
    file1="$1"
    echo "Processing files - download: $file1"

    download_url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$file1/$file1"
    echo "Download URL: $download_url"
    
    wget "$download_url"

    if fastq-dump --split-files --origfmt --gzip "$file1"; then
        if [ $? -ne 0 ]; then
            echo "Error processing files - download: $file1"
        fi
    fi
}

export -f process_files

# Run the process_files function in parallel for each accession number
parallel -j 20 process_files ::: "${files1[@]}"


##################
# Deduplication
##################

mapfile -t files1 < acc_list.txt
printf '%s\n' "${files1[@]}"

process_files() {
    local file1="$1"
    local input_file="${file1}_1.fastq.gz"
    local input_file2="${file1}_2.fastq.gz"
    local output_file="${file1}_1_de.fastq.gz"
    local output_file2="${file1}_2_de.fastq.gz"

    echo "Processing files: $input_file" >> deduplicate.log

    #/work1/home/yenmr/bin/
    if rmDupPE.pl "$input_file" "$input_file2" "$output_file" "$output_file2" >> deduplicate.log; then
        if [ $? -ne 0 ]; then
            echo "Error processing files: $input_file"
        fi
    fi
}

export -f process_files

# Adjust the parallel job count based on your system's capabilities
parallel -j 10 process_files ::: "${files1[@]}"

##################
# moving files
##################

mkdir 1.Raw_Fastq
mv *_de.fastq.gz 1.Raw_Fastq

