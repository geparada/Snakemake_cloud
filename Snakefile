# Snakefile
configfile: "config.yaml"

def get_samples(condition):
    """Return the samples for a given condition."""
    return [sample for sample, cond in config["samples"].items() if cond == condition]

def get_conditions(comparison):
    """Return the conditions for a given comparison."""
    conditions = comparison.split("_vs_")
    conditionA, conditionB = conditions[0], conditions[1]
    return conditionA, conditionB

def deseq2_input_files(wildcards):
    """Return the input files for the deseq2 rule based on the comparison."""
    conditionA, conditionB = get_conditions(wildcards.comparison)
    samplesA = get_samples(conditionA)
    samplesB = get_samples(conditionB)
    all_samples = samplesA + samplesB
    return expand("salmon_quant/{sample}/quant.sf", sample=all_samples)

rule all:
    input:
        expand("results/{comparison}/transcript_differential_expression.tsv", comparison=config["comparisons"]),
        expand("results/{comparison}/gene_differential_expression.tsv", comparison=config["comparisons"])


rule get_fastq_pe:
    output:
        "data/{accession}_1.fastq", "data/{accession}_2.fastq"
    log:
        "logs/pe/{accession}.log"
    params:
        extra="--skip-technical"
    wrapper:
        "v3.13.8/bio/sra-tools/fasterq-dump"

rule salmon_index:
    input:
        sequences=config["transcriptome_fasta"]
    output:
        directory("salmon")
        #"salmon/complete_ref_lens.bin"
        # multiext(
        #     "salmon/",
        #     "complete_ref_lens.bin",
        #     "ctable.bin",
        #     "ctg_offsets.bin",
        #     "duplicate_clusters.tsv",
        #     "info.json",
        #     "mphf.bin",
        #     "pos.bin",
        #     "pre_indexing.log",
        #     "rank.bin",
        #     "refAccumLengths.bin",
        #     "ref_indexing.log",
        #     "reflengths.bin",
        #     "refseq.bin",
        #     "seq.bin",
        #     "versionInfo.json",
        # ),
    params:
        # optional parameters
        extra="--gencode"
    log:
        "logs/salmon_index.log"
    wrapper:
        "v3.13.8/bio/salmon/index"

rule salmon_quant:
    input:
        index="salmon",
        r1="data/{sample}_1.fastq",
        r2="data/{sample}_2.fastq"
    output:
        quant="salmon_quant/{sample}/quant.sf",
        lib="salmon_quant/{sample}/lib_format_counts.json"
    params:
        # optional parameters
        libtype="A",
        extra="",
    wrapper:
        "v3.13.8/bio/salmon/quant"


rule generate_tx2gene:
    input:
        gtf=config["transcriptome_gtf"]
    output:
        tx2gene="salmon_quant/tx2gene.txt"
    shell:
        """
        gzip -dc {input.gtf} | awk -F'\t' '$3 == "transcript" {{ 
            split($9, a, ";"); 
            for(i in a) {{ 
                if(a[i] ~ /gene_id/) {{ gene_id = a[i] }}; 
                if(a[i] ~ /transcript_id/) {{ transcript_id = a[i] }}; 
            }} 
            gsub(/gene_id |"|\s+/, "", gene_id); 
            gsub(/transcript_id |"|\s+/, "", transcript_id); 
            print transcript_id "," gene_id 
        }}' > {output.tx2gene}
        """


rule deseq2_transcript:
    input:
        quants=deseq2_input_files
    output:
        "results/{comparison}/transcript_differential_expression.tsv"
    params:
        conditionA=lambda wildcards: get_samples(get_conditions(wildcards.comparison)[0]),
        conditionB=lambda wildcards: get_samples(get_conditions(wildcards.comparison)[1])
    log:
        "logs/deseq2/transcript_{comparison}.log"
    conda:
        "env/deseq2_env.yaml"
    script:
        "scripts/run_deseq2_transcript_level.R"


rule deseq2_gene:
    input:
        quants=deseq2_input_files,
        tx2gene="salmon_quant/tx2gene.txt"
    output:
        "results/{comparison}/gene_differential_expression.tsv"
    params:
        conditionA=lambda wildcards: get_samples(get_conditions(wildcards.comparison)[0]),
        conditionB=lambda wildcards: get_samples(get_conditions(wildcards.comparison)[1])
    log:
        "logs/deseq2/gene_{comparison}.log"
    conda:
        "env/deseq2_env.yaml"
    script:
        "scripts/run_deseq2_gene_level.R"