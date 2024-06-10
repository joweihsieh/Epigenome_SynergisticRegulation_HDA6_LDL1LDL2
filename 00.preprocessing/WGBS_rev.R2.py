import time
import gzip

def process_fastq(input_file, output_file):
    tStart = time.time()

    with gzip.open(input_file, 'rt') as reader, gzip.open(output_file, 'wt') as writer:
        for index, line in enumerate(reader):
            if index % 4 == 1:
                string = line.strip()
                reversed_string = string[::-1]
                comp_reversed_string = reversed_string.replace("A", "t").replace("T", "a").replace("G", "c").replace("C", "g")
                com_rev_str = comp_reversed_string.upper()
                writer.write(com_rev_str + "\n")
            elif index % 4 == 2:
                writer.write(line.strip() + "\n")
            elif index % 4 == 3:
                qual = line.strip()
                reversed_qual = qual[::-1]
                writer.write(reversed_qual + "\n")
            elif index % 4 == 0:
                writer.write(line.strip() + "\n")

    tEnd = time.time()
    print(f"Processed {input_file} in {tEnd - tStart:.2f} seconds.")

# Example usage with a list of input and output files
input_files = [
    'At-C1_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz',
    'At-C2_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz',
    'At-121_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz',
    'At-122_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz',
    'At-a1_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz',
    'At-a2_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz',
    'At-a121_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz',
    'At-a122_HGL7NCCX2_L8_R2_rmDup.paired.fastq.gz'
]

output_files = [
    'At-C1_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz',
    'At-C2_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz',
    'At-121_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz',
    'At-122_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz',
    'At-a1_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz',
    'At-a2_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz',
    'At-a121_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz',
    'At-a122_HGL7NCCX2_L8_R2_rmDup.paired.rev.fastq.gz'
]

for input_file, output_file in zip(input_files, output_files):
    process_fastq(input_file, output_file)
