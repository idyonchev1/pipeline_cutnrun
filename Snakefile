import pandas as pd
from snakemake.utils import R

# Read the samples file
samples_df = pd.read_csv("samples.txt")
samples = samples_df['Sample'].tolist()

# Create export_name_dict
export_name_dict = {}
sample_from_export = {}

for _, row in samples_df.iterrows():
    sample = row['Sample']
    export_name = f"{row['line']}_{row['condition']}_rep{row['replicate']}"
    export_name_dict[sample] = export_name
    sample_from_export[export_name] = sample

# Function to get non-IgG samples
def get_non_igg_samples():
    return [sample for sample, export_name in export_name_dict.items() if 'IgG' not in export_name]

# Function to get IgG control for a given sample
def get_igg_control(wildcards):
    sample = wildcards.sample.split('_')[0]  # Extract the original sample name
    sample_info = samples_df[samples_df['Sample'] == sample].iloc[0]
    igg_sample = samples_df[(samples_df['condition'] == 'IgG') & 
                            (samples_df['line'] == sample_info['line']) & 
                            (samples_df['replicate'] == sample_info['replicate'])]['Sample'].iloc[0]
    return f"{igg_sample}.bam"

def get_norm_factor(wildcards):
    sample = wildcards.sample.split('_')[0]  # Extract the original sample name
    with open("normalized_counts.tsv") as f:
        next(f)  # Skip the header line
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == sample:
                norm_factor = fields[6]
                print(f"Sample: {sample}, Normalization factor: {norm_factor}")
                return norm_factor
    raise ValueError(f"Sample {sample} not found in normalized_counts.tsv")

rule all:
    input:
        expand("greensheet_norm/{export_name}.bw", export_name=export_name_dict.values()),
        "results.npz",
        "output.tab",
        "normalized_counts.tsv",
        expand("macs2/{sample}_peaks.narrowPeak", 
               sample=[f"{sample}_{export_name}" for sample, export_name in export_name_dict.items() if 'IgG' not in export_name])

#To prevent peakcalling on the drosophila genome we remove those reads, to also prevent FDR dilution on hg38 peaks.
rule filter_bam:
    input:
        bam="{sample}.bam",
        bai="{sample}.bam.bai"
    output:
        "filtered_bam/{sample}.filtered.bam"
    threads: 4
    resources:
        mem_mb=8000,
        time="24:00:00"
    conda:
        "cgat-apps"
    shell:
        """
        samtools view -h {input.bam} | \
        awk '{{if($0 ~ /^@/) {{if($0 !~ /@SQ.*SN:.*_/ && $0 !~ /@SQ.*SN:chrM/) print $0}} else {{if($3 !~ /_/ && $3 != "chrM") print $0}}}}' | \
        samtools view -bS - > {output}
        samtools index {output}
        """

rule get_counts:
    input:
        bams=expand("filtered_bam/{sample}.bam", sample=samples)
    output:
        npz="results.npz",
        counts="output.tab"
    threads: 4
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
        bam="filtered_bam/{sample}.filtered.bam",
        bai="filtered_bam/{sample}.bam.bai",
        normalized="normalized_counts.tsv"
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
        lambda wildcards: f"greensheet_norm/{sample_from_export[wildcards.export_name]}.bedgraph"
    output:
        "greensheet_norm/{export_name}.bw"
    threads: 1
    resources:
        mem_mb=4000,
        time="24:00:00"
    shell:
        "bedGraphToBigWig {input} GRCh38+dm6.chrom.sizes {output} 2>>{output}.log"

rule bam_to_bed:
    input:
        bam="filtered_bam/{sample}.filtered.bam"
    output:
        bed="bed/{sample}.bed"
    threads: 1
    resources:
        mem_mb=8000,
        time="24:00:00"
    conda:
        "cgat-apps"
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output.bed}
        """

rule macs2_peakcalling:
    input:
        treatment = lambda wildcards: f"bed/{wildcards.sample.split('_')[0]}.bed",
        control = lambda wildcards: f"bed/{get_igg_control(wildcards).replace('.bam', '')}.bed",
        normalized_counts = "normalized_counts.tsv"
    output:
        "macs2/{sample}_peaks.narrowPeak"
    resources:
        mem_mb=10000,
        time="24:00:00"
    threads: 1
    params:
        name = "{sample}",
        genome = "hs",
        ratio = get_norm_factor
    log:
        "logs/macs2/{sample}.log"
    conda:
        "cgat-apps"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} \
            -n {params.name} \
            --ratio {params.ratio} \
            --format BED \
            -g {params.genome} \
            --keep-dup all \
            --nomodel \
            -q 0.01 \
            2> {log}
        mv {params.name}_peaks.narrowPeak {output}
        mv {params.name}_peaks.xls macs2
        mv {params.name}_summits.bed macs2
        """
