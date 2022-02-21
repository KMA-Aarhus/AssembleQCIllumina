# snakemake --profile configs/slurm


print("/*")
__author__ = "Tine Ebsen" # Please add your name here if you make changes.
__version__ = "0.3"

from pathlib import Path
import os
import pandas as pd
import sys





print("        /\\                          | |   | |     / __ \\ / ____| KMA,AUH")
print("       /  \\   ___ ___  ___ _ __ ___ | |__ | | ___| |  | | |      ")
print("      / /\\ \\ / __/ __|/ _ \\ '_ ` _ \\| '_ \\| |/ _ \\ |  | | |     ")
print("     / ____ \\__ \\__ \\  __/ | | | | | |_) | |  __/ |__| | |____  ")
print("    /_/    \\_\\___/___/\\___|_| |_| |_|_.__/|_|\\___|\\___\\_\\_____| ")
print(" ")


# Field variables, which might make sense to set from a config file.
#clean_uploads_dir = "../BACKUP/nanopore_sarscov2/pappenheim_clean"
#clean_uploads_dir = "../BACKUP/nanopore_sarscov2/pappenheim_clean/testdir"

out_base = config["out_base"]
sample_reads = config["sample_reads"]

# Set the directory where clean uploads are deposited from the pappenheim workstation pipeline:
in_base = Path(config["raw_uploads_dir"])
kraken2_db = config["kraken2_db"]
##################################
# Parse clean upload directories #
##################################


print("Parsing input directories from", in_base, file = sys.stderr)
files = sorted([str(f) for f in in_base.iterdir() if not f.is_dir()])
print(files)
samples = [s.split("_R1", 1)[0] for s in files if "R1" in s]
R1 = [s for s in files if "R1" in s]
R2 = [s for s in files if "R2" in s]


#df = pd.DataFrame(files, columns = ["filename"])
df = pd.DataFrame(list(zip([s.replace(str(in_base)+"/","") for s in samples], R1, R2)), columns = ["sample_id", "R1", "R2"])
print(df)

rule all:
    input:
        expand(["{out_base}/{sample_id}/trimmed/{sample_id}_val_1.fq.gz", \
                "{out_base}/{sample_id}/trimmed/{sample_id}_val_2.fq.gz", \
                "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz", \
                "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz", \
                "{out_base}/{sample_id}/{sample_id}_kraken2_reads_report.txt", \
                "{out_base}/{sample_id}/{sample_id}_assembly-stats.tab", \
                "{out_base}/{sample_id}/{sample_id}_consensus.fasta", \
                "{out_base}/{sample_id}/{sample_id}.gff", \
                "{out_base}/multiqc_report.html" \
                ], \
                
               out_base = out_base, sample_id = df["sample_id"])





####################
# Setup for data analysis #
####################

if config['trimmed'] == "y":
        # Trim adapters
    rule trim_adapt:
        input: 
            R1 = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["R1"].values[0],
            R2 = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["R2"].values[0]
        output: 
            R1 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_1.fq.gz",
            R2 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_2.fq.gz"
        conda: "configs/conda.yaml"
        threads: 4
        shell: """
        mkdir -p {out_base}/{wildcards.sample_id}/trimmed

        cp {input.R1} {output.R1} 
        cp {input.R2} {output.R2} 
            
        """

else:
    # Trim adapters
    rule trim_adapt:
        input: 
            R1 = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["R1"].values[0],
            R2 = lambda wildcards: df[df["sample_id"]==wildcards.sample_id]["R2"].values[0]
        output: 
            R1 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_1.fq.gz",
            R2 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_2.fq.gz"
        conda: "configs/conda.yaml"
        threads: 4
        shell: """
        mkdir -p {out_base}/{wildcards.sample_id}/trimmed

        trim_galore --paired --gzip --cores 4 --basename {wildcards.sample_id} --fastqc -o {out_base}/{wildcards.sample_id}/trimmed  {input.R1} {input.R2} --length 100 --quality 25
            
        """
## The data I have seen thus far already have UMIs as part of their name:
# @A00606:487:H75CJDSX3:1:1622:10529:12493:AACCACACA 1:N:0:GATCCATG+CAACTCCA where "AACCACACA" is the UMI

# Downsample to an easier-to-handle number of reads
# TODO: Sample size should be replaced with custom values ####
# If samples contains less than the specified value, the output will be identical to the input
rule downsample:
    input: 
        R1 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_1.fq.gz",
        R2 = "{out_base}/{sample_id}/trimmed/{sample_id}_val_2.fq.gz"
    output: 
        R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
        R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz"
    conda: "configs/conda.yaml"
    threads: 1
    shell: """

    seqtk sample -s100 {input.R1} {sample_reads} | gzip -cvf > {output.R1}
    seqtk sample -s100 {input.R2} {sample_reads} | gzip -cvf > {output.R2}
        
    """

#Kraken2
# TODO: Should kraken use the full read set or the downsampled? ####
rule kraken2:
    input:         
        R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
        R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz"
    output: 
        "{out_base}/{sample_id}/{sample_id}_kraken2_reads_report.txt"
    threads: 8
    conda: "configs/conda.yaml"
    shell: """

        kraken2 --db {kraken2_db} --report {output} --threads 8 --paired {input.R1} {input.R2}
        
    """
if config['option'] == 'UMI':

    # Assembly using unicycler
    rule assemble:
        input: 
            R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
            R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz"
        output: 
            contigs = "{out_base}/{sample_id}/consensus/{sample_id}_contigs.fasta",
            assembly_stats = "{out_base}/{sample_id}/{sample_id}_assembly-stats.tab"
        conda: "configs/conda.yaml"
        threads: 8
        shell: """

            mkdir -p {out_base}/{wildcards.sample_id}/assembly

            unicycler --min_fasta_length 500 -1 {input.R1} -2 {input.R2} -o {out_base}/{wildcards.sample_id}/assembly --no_pilon --threads 8
            cp {out_base}/{wildcards.sample_id}/assembly/assembly.fasta {output.contigs}
            assembly-stats -t {output.contigs} > {output.assembly_stats}        
            
            """

    # Map reads to assembly to utilise UMIs

    rule bwa_map:
        input:
            R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
            R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz",
            contigs = "{out_base}/{sample_id}/consensus/{sample_id}_contigs.fasta"
        output:
            "{out_base}/{sample_id}/mapped_reads/{sample_id}.bam"
        conda: "configs/conda.yaml"
        threads: 8
        shell: """
            bwa index {input.contigs}
            bwa mem {input.contigs} {input.R1} {input.R2} -t 8 | samtools sort > {output}
            samtools index {output}

            """

    ## The data I have seen thus far already have UMIs as part of their name:
    # @A00606:487:H75CJDSX3:1:1622:10529:12493:AACCACACA 1:N:0:GATCCATG+CAACTCCA where "AACCACACA" is the UMI        
    rule umi_tools:
        input:
            "{out_base}/{sample_id}/mapped_reads/{sample_id}.bam"
        output:
            "{out_base}/{sample_id}/mapped_reads/{sample_id}_deduplicated.bam"
        conda: "configs/conda.yaml"
        shell: """
            umi_tools dedup -I {input} --output-stats={out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}_output_stats --umi-separator=':' -S {output}
            samtools index {output}

            """

    rule consensus:
        input:
            mapping = "{out_base}/{sample_id}/mapped_reads/{sample_id}_deduplicated.bam",
            contigs = "{out_base}/{sample_id}/consensus/{sample_id}_contigs.fasta"
        output:
            "{out_base}/{sample_id}/{sample_id}_consensus.fasta"
        conda: "configs/conda.yaml"
        shell: """

            bcftools mpileup --fasta-ref {input.contigs} {input.mapping} | bcftools call -m -o {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf
            bgzip -f {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf
            tabix {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf.gz
            bcftools consensus --fasta-ref {input.contigs} {out_base}/{wildcards.sample_id}/mapped_reads/{wildcards.sample_id}.vcf.gz -o {output}

            """

else:
    # Assembly using unicycler
    rule assemble:
        input: 
            R1 = "{out_base}/{sample_id}/sampled/{sample_id}_R1_sampled.fq.gz",
            R2 = "{out_base}/{sample_id}/sampled/{sample_id}_R2_sampled.fq.gz"
        output: 
            contigs = "{out_base}/{sample_id}/{sample_id}_consensus.fasta",
            assembly_stats = "{out_base}/{sample_id}/{sample_id}_assembly-stats.tab"
        conda: "configs/conda.yaml"
        threads: 8
        shell: """

            mkdir -p {out_base}/{wildcards.sample_id}/assembly

            unicycler --min_fasta_length 500 -1 {input.R1} -2 {input.R2} -o {out_base}/{wildcards.sample_id}/assembly --threads 8
            cp {out_base}/{wildcards.sample_id}/assembly/assembly.fasta {output.contigs}
            assembly-stats -t {output.contigs} > {output.assembly_stats}        
            
            """

rule annotate_genes:
    input:
        "{out_base}/{sample_id}/{sample_id}_consensus.fasta"
    output:
        "{out_base}/{sample_id}/{sample_id}.gff"
    conda: "configs/prokka.yaml"
    threads: 8
    shell: """
        mkdir -p {out_base}/{wildcards.sample_id}/prokka
        prokka --outdir {out_base}/{wildcards.sample_id}/prokka --cpu 8 --force --prefix {wildcards.sample_id} {input}
        cp {out_base}/{wildcards.sample_id}/prokka/{wildcards.sample_id}.gff {out_base}/{wildcards.sample_id}/{wildcards.sample_id}.gff 

        """

rule multiqc:
    input:
        expand("{out_base}/{sample_id}/{sample_id}.gff", out_base = out_base, sample_id = df["sample_id"])
    output:
        "{out_base}/multiqc_report.html"
    conda: "configs/conda.yaml"
    threads: 1
    shell: """
        multiqc -d {out_base} -o {out_base}

        """

