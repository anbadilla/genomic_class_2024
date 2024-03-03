# Importing necessary modules
import os

# Define input and output directories
input_dir = "/afs/crc.nd.edu/user/a/abadilla/genomics_project/raw_data"
output_dir = "/afs/crc.nd.edu/user/a/abadilla/genomics_project/results"
kraken_db = "/afs/crc.nd.edu/user/a/abadilla/rob_sequencing/kraken2_scripts/kraken2-master/db/silva138"

# Define a list of sample names
samples = os.listdir(input_dir)

# Rule to run FastQC on each sample
rule fastqc:
    input:
        fastq1 = "{input_dir}/{sample}_R1.fastq.gz",
        fastq2 = "{input_dir}/{sample}_R2.fastq.gz"
    output:
        html_report = "{output_dir}/{sample}_fastqc.html",
        zip_file = "{output_dir}/{sample}_fastqc.zip"
    shell:
        "fastqc {input.fastq1} {input.fastq2} --outdir {output_dir}"

# Rule to run Trimmomatic to trim adapters and low-quality bases
rule trimmomatic:
    input:
        fastq1 = "{input_dir}/{sample}_R1.fastq.gz",
        fastq2 = "{input_dir}/{sample}_R2.fastq.gz"
    output:
        trimmed_fastq1 = "{output_dir}/{sample}_R1.trimmed.fastq.gz",
        trimmed_fastq2 = "{output_dir}/{sample}_R2.trimmed.fastq.gz"
    shell:
        "trimmomatic PE -threads 4 {input.fastq1} {input.fastq2} "
        "{output.trimmed_fastq1} {output_dir}/{sample}_R1_unpaired.fastq.gz "
        "{output.trimmed_fastq2} {output_dir}/{sample}_R2_unpaired.fastq.gz "
        "ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Rule to run Kraken2 and Bracken
rule kraken_bracken:
    input:
        fastq1 = "{input_dir}/{sample}_R1.fastq.gz",
        fastq2 = "{input_dir}/{sample}_R2.fastq.gz"
    output:
        kraken_report = "{output_dir}/{sample}.kraken2",
        bracken_report = "{output_dir}/{sample}.bracken"
    shell:
        "kraken2 --db {kraken_db} --threads 1 --minimum-hit-groups 3 --report-minimizer-data --report {output.kraken_report} --paired {input.fastq1} {input.fastq2} > {output.kraken_report}; "
        "bracken -d {kraken_db} -i {output.kraken_report} -o {output.bracken_report} -w {output.bracken_report}.kreport2 -r 250 -l G"

# Rule to define all outputs
rule all:
    input:
        expand("{output_dir}/{sample}.kraken2", output_dir=output_dir, sample=samples),
        expand("{output_dir}/{sample}.bracken", output_dir=output_dir, sample=samples)

