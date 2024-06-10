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

