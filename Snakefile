#shell.prefix("set -euo pipefail;") 

configfile: "config.yaml"

# Folder variables
assembly     = "data/assembly"
raw          = "data/fastq_raw"
trimmed      = "data/fastq_trimmed"
indexes      = "data/indexes" 
quant        = "data/quant"
sleuth       = "data/sleuth"
transdecoder = "data/transdecoder"
filtering    = "data/filtering"
g2t_k2       = "data/g2t_k2"


# Path to programs
kallisto    = config["software"]["kallisto"]
trimmomatic = config["software"]["trimmomatic"]
samtools    = config["software"]["samtools"]
gzip        = config["software"]["gzip"]
samtools    = config["software"]["samtools"]
bcftools    = config["software"]["bcftools"]
bowtie2build=   "bowtie2-build"
bowtie2     = "bowtie2"
samtools0   = "./bin/samtools0"
bcftools0   = "./bin/bcftools0"

rule all:
    input:
        g2t_k2 + "/all.vcf",
        expand(
            quant + "/{sample}/{sample}.bam",
            sample = config["rna_pe"]
        )
        

rule clean:
    shell:
        """
        rm -rf data/fastq_trimmed
        rm -rf data/indexes
        rm -rf data/quant
        rm -rf data/sleuth
        rm -rf data/filtering
        rm -rf data/g2t_k2
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
        "benchmarks/kallisto/build_kallisto_index.json"
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
        pseudobam     = temp(quant + "/{sample}/{sample}_unsorted.bam")
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
            -s - -                      \
        > {output.pseudobam} )          \
        2> {log}
        """



rule kallisto_sort_pseudobam:
    input:
        pseudobam = quant + "/{sample}/{sample}_unsorted.bam"
    output:
        pseudobam = quant + "/{sample}/{sample}.bam"
    threads:
        24
    log:
        "logs/quant/sort_pseudobam_{sample}.log"
    benchmark:
        "benchmarks/quant/sort_pseudobam_{sample}.json"
    shell:
        """
        {samtools} sort             \
            -@  {threads}           \
            -l  9                   \
            -o  {output.pseudobam}  \
            -T  $(mktemp -d)        \
            -O  bam                 \
            {input.pseudobam}       \
        2>  {log}
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
        "benchmarks/sleuth/get_normalized_tpm_table.json"
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
        "benchmarks/filtering/get_expressed_isoforms_ids.json"
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
        "benchmarks/filtering/get_coding_isoforms_ids.json"
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
        "benchmarks/filtering/get_expressed_and_coding_ids.json"
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
        "benchmarks/filtering/get_expressed_and_coding_fasta.log"
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
        "benchmarks/filtering/compute_monotigs_id.json"
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
        "benchmarks/filtering/get_monotigs_fasta.json"
    shell:
        """
        samtools faidx {input.filtered_assembly}
        
        cat {input.monotig_ids}                         |
        xargs samtools faidx {input.filtered_assembly}  \
        >  {output.monotig_assembly}                    \
        2> {log}
        """


rule create_monotigs_index:
    input:
        fasta= filtering + "/id_expressed_and_coding.fasta"
    output:
        index= filtering + "/id_expressed_and_coding.fasta.fai"
    threads:
        1
    log:
        "logs/filtering/create_monotigs_index.log"
    benchmark:
        "benchmarks/filtering/create_monotigs_index.log"



#############
## G2T k=2 ##
#############

rule g2t_k2_make_bowtie2_index_assembly_exp_cod_mono:
    input:
        monotig_assembly= filtering + "/assembly_exp_cod_mono.fasta"  
    output:
        aux= expand(indexes + "/assembly_exp_cod_mono.{extension}",
            extension= "1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2".split()),
        name= touch(indexes + "/assembly_exp_cod_mono")
    threads:
        1
    log:
        "logs/g2t_k2/g2t_k2_make_bowtie2_index_assembly_exp_cod_mono.log"
    benchmark:
        "benchmarks/g2t_k2/g2t_k2_make_bowtie2_index_assembly_exp_cod_mono"
    shell:
        """
        {bowtie2build}                  \
            {input.monotig_assembly}    \
            {output.name}               \
        >   {log} 2>&1
        """
    



rule g2t_k2_mapping:
    """
    Perform the g2t mapping with bowtie2 and the 
    sam | bam | rmdup | sort conversions
    """
    input:
        index=  indexes + "/assembly_exp_cod_mono",
        f=      trimmed + "/{sample}_1.fastq.gz",
        r=      trimmed + "/{sample}_2.fastq.gz",
        u=      trimmed + "/{sample}_u.fastq.gz"
    output:
        bam=    g2t_k2 + "/{sample}.bam"
    log:
        "logs/g2t_k2/g2t_k2_mapping_{sample}.log"
    benchmark:
        "benchmarks/g2t_k2/g2t_k2_mapping_{sample}.json"
    params:
        id =      lambda wildcards: wildcards.sample,
        library = lambda wildcards: "LB:truseq_" + wildcards.sample,
        platform= "PL:Illumina",
        sample=   lambda wildcards: "SM:" + wildcards.sample
    threads:
        24
    shell:
        """
        ( {bowtie2}                         \
            --threads   {threads}           \
            --phred33                       \
            --quiet                         \
            --sensitive-local               \
            --no-unal                       \
            --no-unal                       \
            --rg-id     {params.id}         \
            --rg        {params.library}    \
            --rg        {params.platform}   \
            --rg        {params.sample}     \
            -k  2                           \
            -x  {input.index}               \
            -1  {input.f}                   \
            -2  {input.r}                   \
            -U  {input.u}                   |
        {samtools} view                     \
            -S  `# Input is SAM`            \
            -h  `# Print SAM header`        \
            -u  `# Uncompressed SAM`        \
            -                               |
        {samtools0} rmdup -s - -            |
        {samtools} sort                     \
            -@  {threads}                   \
            -l  9                           \
            -o  -                           \
            -T  $(mktemp -d)                \
            -O  bam                         \
        > {output.bam} )                    \
        2> {log}
        """
        




rule g2t_k2_index_bam:
    input:
        bam = g2t_k2 + "/{sample}.bam"
    output:
        bai = g2t_k2 + "/{sample}.bam.bai"
    threads:
        1
    log:
        "logs/g2t_k2/g2t_k2_index_bam_{sample}.log"
    benchmark:
        "benchmarks/g2t_k2/g2t_k2_index_bam_{sample}.json"
    shell:
        """
        {samtools} index {input.bam} 2> {log}
        """



rule g2t_k2_mpileup:
    input:
        assembly= filtering + "/assembly_exp_cod_mono.fasta",
        bams= expand(
            g2t_k2 + "/{sample}.bam",
            sample = config["dna_pe"]
        ),
        indexes= expand(
            g2t_k2 + "/{sample}.bam.bai",
            sample = config["dna_pe"]
        )
    output:
        bcf= g2t_k2 + "/all.bcf"
    threads:
        1
    log:
        "logs/g2t_k2/g2t_k2_mpileup.log"
    benchmark:
        "benchmarks/g2t_k2/g2t_k2_mpileup.json"
    shell:
        """
        {samtools0} mpileup     \
        -ABDguI                 \
        -f  {input.assembly}    \
        -C  50                  \
        {input.bams}            \
        >   {output.bcf}        \
        2>  {log}
        """



rule g2t_k2_make_samples_txt:
    output:
        tsv= g2t_k2 + "/samples.txt"
    threads:
        1
    log:
        "logs/g2t_k2/g2t_k2_make_samples_txt.log"
    benchmark:
        "benchmarks/g2t_k2/g2t_k2_make_samples_txt.json"
    run:
        SAMPLES = config["dna_pe"]
        with open(output.tsv, "w") as f_out:
            for sample in SAMPLES:
                f_out.write(sample + "\t" + "2" + "\n")



rule g2t_k2_bcf_to_vcf:
    """
    ####################DIRTY
    """
    input:
        bcf=     g2t_k2 + "/all.bcf",
        samples= g2t_k2 + "/samples.txt"
    output:
        vcf=     g2t_k2 + "/all.vcf"
    threads:
        1
    log:
        "logs/g2t_k2/g2t_k2_bcf_to_vcf.log"
    benchmark:
        "benchmarks/g2t_k2/g2t_k2_bcf_to_vcf.json"
    shell:
        """
        {bcftools0} view -LNcegIv   \
            -s {input.samples}      \
            {input.bcf}             \
        >   {output.vcf}            \
        2>  {log}
        """
