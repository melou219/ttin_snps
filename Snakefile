shell.prefix("set -euo pipefail;") 

configfile: "config.yaml"

# Folder variables
raw     = "data/fastq_raw"
trimmed = "data/fastq_trimmed"
indexes = "data/indexes" 
quant   = "data/kallisto_quant"



# Path to programs
kallisto    = config["software"]["kallisto"]
trimmomatic = config["software"]["trimmomatic"]
samtools    = config["software"]["samtools"]
gzip        = config["software"]["gzip"]
samtools    = config["software"]["samtools"]
bcftools    = config["software"]["bcftools"]
bowtie2build=   "bowtie2-build"


rule all:
    input:
        expand(trimmed + "/{sample}_u.fastq.gz",
            sample= config["dna_pe"]),
        expand(trimmed + "/{sample}_u.fastq.gz",
            sample= config["rna_pe"]),
        config["assembly"] + ".fai",
        indexes + "/assembly_kallisto.idx",
        expand(quant + "/{sample}/abundance.tsv",
            sample= config["rna_pe"])



rule clean:
    shell:
        """
        rm -rf data/fastq_g_trimmed
        rm -rf data/fastq_t_trimmed
        rm -rf data/indexes
        rm -rf data/kallisto_quant
        rm -rf benchmarks
        rm -rf logs
        """



rule index_assembly_samtools:
    input:
        assembly=   config["assembly"]
    output:
        index=      config["assembly"] + ".fai"
    log:
        "logs/assembly/index.log"
    benchmark:
        "benchmarks/index_assembly.json"
    threads:
        1
    shell:
        """
        {samtools} faidx {input.assembly}   \
        >   {log}
        """



# Trim reads
rule trimmomatic_rna:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and remove
    low quality regions and reads.
    """
    input:
        forward    = lambda wildcards: config["rna_pe"][wildcards.sample]["forward"],
        reverse    = lambda wildcards: config["rna_pe"][wildcards.sample]["reverse"]
    output:
        forward    = protected("data/fastq_trimmed/{sample}_1.fastq.gz"),
        reverse    = protected("data/fastq_trimmed/{sample}_2.fastq.gz"),
        unpaired   = protected("data/fastq_trimmed/{sample}_u.fastq.gz")
    params:
        unpaired_1 = trimmed + "/{sample}_3.fastq.gz",
        unpaired_2 = trimmed + "/{sample}_4.fastq.gz",
        adaptor    = lambda wildcards: config["rna_pe"][wildcards.sample]["adaptor"],
        trimming   = config["trimmomatic_params"],
        phred      = lambda wildcards: config["rna_pe"][wildcards.sample]["phred"]
    benchmark:
        "benchmarks/trimmomatic/{sample}.json"
    log:
        "logs/trimmomatic/{sample}.log" 
    threads:
        4 # It doesn't perform well above this value
    shell:
        """
        {trimmomatic} PE                            \
            -threads {threads}                      \
            -{params.phred}                         \
            {input.forward}                         \
            {input.reverse}                         \
            {output.forward}                        \
            {params.unpaired_1}                     \
            {output.reverse}                        \
            {params.unpaired_2}                     \
            ILLUMINACLIP:{params.adaptor}:2:30:10   \
            {params.trimming}                       \
            2> {log}
            
        {gzip} -dc {params.unpaired_1} {params.unpaired_2}  |
        {gzip} -9 > {output.unpaired}
        
        rm {params.unpaired_1} {params.unpaired_2}
        """


rule trimmomatic_dna:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and remove
    low quality regions and reads.
    """
    input:
        forward    = lambda wildcards: config["dna_pe"][wildcards.sample]["forward"],
        reverse    = lambda wildcards: config["dna_pe"][wildcards.sample]["reverse"]
    output:
        forward    = protected(trimmed + "/{sample}_1.fastq.gz"),
        reverse    = protected(trimmed + "/{sample}_2.fastq.gz"),
        unpaired   = protected(trimmed + "/{sample}_u.fastq.gz")
    params:
        unpaired_1 = trimmed + "/{sample}_3.fastq.gz",
        unpaired_2 = trimmed + "/{sample}_4.fastq.gz",
        adaptor    = lambda wildcards: config["dna_pe"][wildcards.sample]["adaptor"],
        trimming   = config["trimmomatic_params"],
        phred      = lambda wildcards: config["dna_pe"][wildcards.sample]["phred"]
    benchmark:
        "benchmarks/trimmomatic/{sample}.json"
    log:
        "logs/trimmomatic/{sample}.log" 
    threads:
        4 # It doesn't perform well above this value
    shell:
        """
        {trimmomatic} PE                            \
            -threads {threads}                      \
            -{params.phred}                         \
            {input.forward}                         \
            {input.reverse}                         \
            {output.forward}                        \
            {params.unpaired_1}                     \
            {output.reverse}                        \
            {params.unpaired_2}                     \
            ILLUMINACLIP:{params.adaptor}:2:30:10   \
            {params.trimming}                       \
            2> {log}
            
        {gzip} -dc {params.unpaired_1} {params.unpaired_2}  |
        {gzip} -9 > {output.unpaired}
        
        rm {params.unpaired_1} {params.unpaired_2}
        """



rule build_kallisto_index:
    """
    Self-explanatory
    """
    input:
        assembly = config["assembly"]
    output:
        index    = indexes + "/assembly_kallisto.idx"
    benchmark:
        "benchmark/kallisto/build_kallisto_index.json"
    log:
        "logs/kallisto/build_kallisto_index.json"
    threads:
        1
    shell:
        """
        {kallisto} index            \
            --index={output.index}  \
            {input.assembly}        \
        2> {log}
        """



rule kallisto_quant_pe:
    input:
        forward = trimmed + "/{sample}_1.fastq.gz",
        reverse = trimmed + "/{sample}_2.fastq.gz",
        index   = indexes + "/assembly_kallisto.idx"
    output:
        abundance_tsv = quant + "/{sample}/abundance.tsv",
        abundance_h5  = quant + "/{sample}/abundance.h5",
        pseudobam     = quant + "/{sample}/{sample}.bam"
    params:
        outdir  = quant + "/{sample}",
        params  = config["kallisto_params"]
    threads:
        24
    log:
        "logs/quant/{sample}.log"
    benchmark:
        "benchmarks/quant/{sample}.json"
    shell:
        """
        ({kallisto} quant               \
            --index={input.index}       \
            --output={params.outdir}    \
            --pseudobam                 \
            --threads={threads}         \
            {params.params}             \
            {input.forward}             \
            {input.reverse}             |
        {samtools} view                 \
            -S  `# Input is SAM`        \
            -h  `# Print SAM header`    \
            -u  `# Uncompressed SAM`    \
            -                           |
        {samtools} rmdup                \
            -s - -                      |
        {samtools} sort                 \
            -@  {threads}               \
            -l  9                       \
            -o  -                       \
            -T  $(mktemp)               \
            -O  bam                     \
        > {output.pseudobam} )          \
        2> {log}
        """



#############
## G2T k=2 ##
#############

rule make_bowtie2_index:
    input:
        assembly=   config["assembly"]
    output:
        name=       touch("data/index/assembly_bwt")
    log:
        "logs/bowtie2/make_bowtie2_index.log"
    threads:
        1
    params:
        dir=        "data/index"
    shell:
        """
        mkdir -p {params.dir}
    
        {bowtie2build}          \
            {input.assembly}    \
            {output}            \
        >   {log}
        """
    



rule g2t_k2_mapping:
    """
    Perform the g2t mapping with bowtie2 and the 
    sam | bam | rmdup | sort conversions
    """
    input:
        index=  "data/index/assembly",
        f=      "data/fastq_g_trimmed/{sample}_1.fastq.gz",
        r=      "data/fastq_g_trimmed/{sample}_2.fastq.gz",
        u=      "data/fastq_g_trimmed/{sample}_U.fastq.gz"
    output:
        bam=    "data/g2t_k2/{sample}.bam"
    log:
        out=    "data/g2t_k2/{sample}.out",
        err=    "data/g2t_k2/{sample}.err"
    params:
        sample= "{sample}",
        dir=    "data/g2t_k2"
    threads:
        24
    shell:
        """
        mkdir -p {params.dir}
        
        {bowtie2}                   \
            --threads   {threads}   \
            --phred33               \
            --quiet                 \
            --sensitive-local       \
            --no-unal               \
            -k  2                   \
            -x  {input.index}       \
            -1  {input.f}           \
            -2  {input.r}           \
            -U  {input.u}           |
        {samtools}  view            \
            -S  `# Input is SAM`    \
            -h  `# Print SAM header`\
            -u  `# Uncompressed SAM`\
            -                       |
        {samtools}  sort            \
            -@  {threads}           \
            -l  0                   \
            -o                      \
            - $(mktemp)             |
        {samtools}  rmdup           \
            -                       \
        {output.bam}                \
        2>  {log.err}

        touch {log.out}
        """

