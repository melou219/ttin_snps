shell.prefix("set -euo pipefail;") 

configfile: "config.yaml"

# Folder variables
assembly     = "data/assembly"
raw          = "data/fastq_raw"
trimmed      = "data/fastq_trimmed"
indexes      = "data/indexes" 
quant        = "data/kallisto_quant"
sleuth       = "data/sleuth"
transdecoder = "data/transdecoder"
filtering    = "data/filtering"

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
            sample= config["dna_pe"]), # Trimmed genmoe
         filtering + "/assembly_exp_cod_mono.fasta"



rule clean:
    shell:
        """
        rm -rf data/fastq_g_trimmed
        rm -rf data/fastq_t_trimmed
        rm -rf data/indexes
        rm -rf data/kallisto_quant
        rm -rf data/sleuth
        rm -rf data/filtering
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



rule create_experimental_design_for_sleuth:
    output:
        table = sleuth + "/experimental_design.txt"
    threads:
        1
    log:
        "logs/sleuth/create_experimental_design_for_sleuth.log"
    benchmark:
        "benchmarks/sleuth/create_experimental_design_for_sleuth.json"
    run:
        design_dict = {accession: config["rna_pe"][accession]["condition"] for accession in config["rna_pe"]}            
        with open(output.table, "w") as table:
            table.write("run_accession\tcondition\n")
            for accession in design_dict:
                table.write(accession + "\t" + design_dict[accession] + "\n")
            


rule get_normalized_tpm_table:
    input:
        design = sleuth + "/experimental_design.txt",
        tsvs   = expand(
            quant + "/{sample}/abundance.tsv",
                sample = config["rna_pe"]
            ),
        h5s    = expand(
            quant + "/{sample}/abundance.h5",
                sample = config["rna_pe"]
            ),
        bams   = expand(
            quant + "/{sample}/{sample}.bam",
                sample = config["rna_pe"]
            )
    output:
        table= sleuth + "/tpm_normalized.tsv",
        rdata= sleuth + "/sleuth_object.RData"
    threads:
        1
    log:
        "logs/sleuth/get_normalized_tpm_table.log"
    benchmark:
        "benchmark/sleuth/get_normalized_tpm_table.json"
    shell:
        """
        Rscript scripts/sleuth_get_normalized_tpm_table.R \
        2> {log}   
        """



rule get_expressed_isoforms_ids:
    input:
        table= sleuth + "/tpm_normalized.tsv"
    output:
        expressed = filtering + "/id_expressed.tsv"
    threads:
        1
    log:
        "logs/filtering/get_expressed_isoforms_ids.log"
    benchmark:
        "benchmark/filtering/get_expressed_isoforms_ids.json"
    shell:
        """
        ( tail -n +2 {input.table}  |
        awk '$4 >0'                 |
        cut -f 1                    | 
        uniq                        \
        > {output.expressed} )      \
        2> {log}
        """



rule get_coding_isoforms_ids:
    input:
        pep    = config["pep"]
    output:
        coding = filtering + "/id_coding.tsv" 
    threads:
        1
    log:
        "logs/fitering/get_coding_isoforms_ids.log"
    benchmark:
        "benchmark/filtering/get_coding_isoforms_ids.json"
    shell:
        """
        python3 scripts/fasta_to_id.py < {input.pep}    |
        python3 scripts/pep_id_to_isoform.py            |
        sort -u                                         \
        >  {output.coding}                              \
        2> {log}
        """


rule get_expressed_and_coding_ids:
    input:
        expressed = filtering + "/id_expressed.tsv",
        coding    = filtering + "/id_coding.tsv",
    output:
        filter    = filtering + "/id_expressed_and_coding.tsv"
    log:
        "logs/filtering/get_expressed_and_coding_ids.log"
    benchmark:
        "benchmark/filtering/get_expressed_and_coding_ids.json"
    threads:
        1
    shell:
        """
        comm -12                        \
            <(sort {input.expressed})   \
            <(sort {input.coding})      \
        > {output.filter}               \
        2>  {log}
        """


rule get_expressed_and_coding_fasta:
    input:
        assembly= config["assembly"],
        index=    config["assembly"] + ".fai",
        ids=      filtering + "/id_expressed_and_coding.tsv"
    output:
        assembly= filtering + "/id_expressed_and_coding.fasta"
    threads:
        1
    log:
        "logs/filtering/get_expressed_and_coding_fasta.log"
    benchmark:
        "benchmark/filtering/get_expressed_and_coding_fasta.log"
    shell:
        """
        cat {input.ids}                         |
        xargs samtools faidx {input.assembly}   \
        >  {output.assembly}                    \
        2> {log}
        """



rule get_monotigs_id:
    input:
        filtered_assembly= filtering + "/id_expressed_and_coding.fasta"
    output:
        monotig_ids= filtering + "/id_monotigs.tsv"
    threads:
        1
    log:
        "logs/filtering/compute_monotigs_id.log"
    benchmark:
        "benchmark/filtering/compute_monotigs_id.json"
    shell:
        """
        python3 scripts/filter_n_isogroups.py   \
            <(python3 scripts/fasta_to_id.py    \
                < {input.filtered_assembly} )   \
            {output.monotig_ids}                \
        > {log}
        """



rule get_monotigs_fasta:
    input:
        filtered_assembly= filtering + "/id_expressed_and_coding.fasta",
        monotig_ids= filtering + "/id_monotigs.tsv"
    output:
        monotig_assembly= filtering + "/assembly_exp_cod_mono.fasta"
    threads:
        1
    log:
        "logs/filtering/get_monotigs_fasta.log"
    benchmark:
        "benchmark/filtering/get_monotigs_fasta.json"
    shell:
        """
        samtools faidx {input.filtered_assembly}
        
        cat {input.monotig_ids}                         |
        xargs samtools faidx {input.filtered_assembly}  \
        >  {output.monotig_assembly}                    \
        2> {log}
        """


#############
## G2T k=2 ##
#############

rule make_bowtie2_index:
    input:
        monotig_assembly= filtering + "/assembly_exp_cod_mono.fasta"  
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

