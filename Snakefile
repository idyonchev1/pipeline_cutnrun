import pandas as pd
from snakemake.utils import R

# Read the samples file
samples_df = pd.read_csv("samples.txt")
samples = samples_df['Sample'].tolist()

# Create export_name_dict
export_name_dict = {}

for _, row in samples_df.iterrows():
    sample = row['Sample']
    export_name = f"{row['line']}_{row['condition']}_rep{row['replicate']}"
    export_name_dict[sample] = export_name

# Create a reverse mapping from export_name to sample
sample_from_export = {v: k for k, v in export_name_dict.items()}

def get_norm_factor(wildcards):
    with open("normalized_counts.tsv") as f:
        next(f)  # Skip the header line
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == wildcards.sample:
                norm_factor = fields[6]
                print(f"Sample: {wildcards.sample}, Normalization factor: {norm_factor}")
                return norm_factor
    raise ValueError(f"Sample {wildcards.sample} not found in normalized_counts.tsv")

# Function to get the export name for a sample
def get_export_name(wildcards):
    return export_name_dict[wildcards.sample]

# Function to get the sample name from export name
def get_sample_from_export(wildcards):
    return sample_from_export[wildcards.export_name]

rule all:
    input:
        expand("greensheet_norm/{export_name}.bw", export_name=export_name_dict.values()),
        "results.npz",
        "output.tab",
        "normalized_counts.tsv"

rule get_counts:
    input:
        bams=expand("{sample}.bam", sample=samples)
    output:
        npz="results.npz",
        counts="output.tab"
    threads: 4  # Adjust as needed
    resources:
        mem_mb=16000,
        time="24:00:00"
    conda:
        "cgat-apps"
    shell:
        "multiBamSummary bins --bamfiles {input.bams} "
        "--BED hg38_greenlist.bed --smartLabels -e --centerReads "
        "-o {output.npz} --outRawCounts {output.counts}"

rule process_counts:
    input:
        counts="output.tab",
        samples="samples.txt"
    output:
        normalized="normalized_counts.tsv"
    threads: 1
    resources:
        mem_mb=4000,
        time="24:00:00"
    conda:
        "r_env"
    script:
        "process_counts.R"

rule export_bedgraph_temp:
    input:
        bam="{sample}.bam",
        bai="{sample}.bam.bai",
        normalized="normalized_counts.tsv"  # Add this to ensure normalization is done
    output:
        temp("greensheet_norm/{sample}.bedgraph")
    threads: 1
    resources:
        mem_mb=8000,
        time="24:00:00"
    params:
        norm_factor=get_norm_factor
    shell:
        "genomeCoverageBed -bg -ibam {input.bam} -scale {params.norm_factor} | sort -k1,1 -k2,2n > {output}"

rule export_bigwig:
    input:
        lambda wildcards: f"greensheet_norm/{get_sample_from_export(wildcards)}.bedgraph"
    output:
        "greensheet_norm/{export_name}.bw"
    threads: 1
    resources:
        mem_mb=4000,
        time="24:00:00"
    shell:
        "bedGraphToBigWig {input} GRCh38+dm6.chrom.sizes {output} 2>>{output}.log"
