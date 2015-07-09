shell.prefix("set -euo pipefail;") 

gsamples = """
    C6KB3ANXX_3_10nf
    C6KB3ANXX_8_14nf
    """.split()

tsamples = """
    C6CGAACXX_8_14
    C6CGAACXX_8_15
    C6CGAACXX_8_16
    C6CGAACXX_8_18
    C70VFACXX_7_19
    C72MEACXX_7_19
    C72U8ACXX_3_19
    C72U8ACXX_7_19
    """.split()

ends= "1 2".split()
endsu= "1 2 U".split()



## Software
trimmomatic=    "java -jar ./src/Trimmomatic-0.33/trimmomatic-0.33.jar"
bowtie2=        "bowtie2"
bowtie2build=   "bowtie2-build"
samtools=       "./bin/samtools" # Harvest of 2011


rule all:
    input:
        "data/assembly/assembly.fasta",
        expand("data/fastq_g_raw/{sample}_{end}.fastq.gz",
            sample= gsamples,
            end=    ends ),
        expand("data/fastq_t_raw/{sample}_{end}.fastq.gz",
            sample= tsamples,
            end=    ends ),
        expand("data/fastq_g_trimmed/{sample}_{end}.fastq.gz",
            sample= gsamples,
            end=    endsu ),
        expand("data/fastq_t_trimmed/{sample}_{end}.fastq.gz",
            sample= tsamples,
            end=    endsu ),
        "data/index/assembly",
        expand("data/g2t_k2/{sample}.bam",
            sample= gsamples)



rule clean:
    shell:
        """
        rm -rf data/assembly
        rm -rf data/fastq_g_raw
        rm -rf data/fastq_t_raw
        rm -rf data/fastq_g_trimmed
        rm -rf data/fastq_t_trimmed
        """


rule link_transcriptome_assembly:
    """
    Make a link from somewhere in the server (the ttin_assembly project) and 
    leave the link in the data/assembly folder.
    """
    input:
        assembly=   "/home/SHARE/ttin_assembly/data/assembly/all_normalized.fasta"
    output:
        assembly=   "data/assembly/assembly.fasta"
    params:
        folder=     "data/assembly"
    threads:
        1
    log:
        # No log
    shell:
        """
        mkdir -p  {params.folder}
        ln -s {input.assembly} {output.assembly}
        """



rule link_fastq_g_raw:
    """
    Make a link of all the genomic raw fastq.gz files that the CNAG gave us and 
    leave them in the data/fastq_g_raw folder.
    """
    input:
        fastqgz=    "/home/SHARE/raw_data/ttin/ESTONBA_03/20150526/FASTQ/{sample}_{end}.fastq.gz"
    output:
        fastqgz=    "data/fastq_g_raw/{sample}_{end}.fastq.gz"
    threads:
        1
    log:
        # No log
    params:
        dir=        "data/fastq_g_raw"
    shell:
        """
        mkdir -p {params.dir}
        ln -s   {input.fastqgz} {output.fastqgz}
        """



rule link_fastq_t_raw:
    """
    Make a link of all the transcriptomic raw fastq.gz files that the CNAG gave
    us and leave them in the data/fastq_t_raw folder.
    """
    input:
        fastqgz=    "/home/SHARE/raw_data/ttin/ESTONBA_02/20150630/FASTQ/{sample}_{end}.fastq.gz"
    output:
        fastqgz=    "data/fastq_t_raw/{sample}_{end}.fastq.gz"
    threads:
        1
    log:
        # No log
    params:
        dir=    "data/fastq_t_raw"
    shell:
        """
        mkdir -p {params.dir}
        ln -s {input.fastqgz} {output.fastqgz}
        """



rule link_fastq_t_trimmed:
    """
    Make a link of all the trimmed4 files that I did when the ttin_assembly 
    project. Save them in the data/fastq_t_trimmed folder.
    """
    input:
        fastqgz=    "/home/SHARE/ttin_assembly/data/fastq_trimmed/{sample}_{end}.trimmed4.fastq.gz"
    output:
        fastqgz=    "data/fastq_t_trimmed/{sample}_{end}.fastq.gz"
    params:
        dir=        "data/fastq_t_trimmed"
    log:
        # No log
    threads:
        1
    shell:
        """
        mkdir -p {params.dir}
        ln -s {input.fastqgz} {output.fastqgz}
        """



rule trimmomatic_g:
    """
    Trim the genome with trimmomatic. Params are the same as transcriptome:
        HEADCROP:13  - Delete first 13 nt (Random hexamer priming Dudoit et al.)        
        ILLUMINACLIP - Search and remove adaptors (Truseq3-PE2)
        MINLEN:31    -  
        AVGQUAL:10   -
        MINLEN:31    -
        TRAILING:19  -
        MINLEN:31    -
        TOPHRED33    -
    """
    input:
        f=          "data/fastq_g_raw/{sample}_1.fastq.gz",
        r=          "data/fastq_g_raw/{sample}_2.fastq.gz"
    output:
        f=          "data/fastq_g_trimmed/{sample}_1.fastq.gz",
        r=          "data/fastq_g_trimmed/{sample}_2.fastq.gz",
        u=          "data/fastq_g_trimmed/{sample}_U.fastq.gz"
    params:
        dir=        "data/fastq_g_trimmed",
        adaptors=   "./src/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa",
        u1=         "data/fastq_g_trimmed/{sample}_u1.fastq.gz",
        u2=         "data/fastq_g_trimmed/{sample}_u2.fastq.gz"
    threads:
        4
    log:
        out=        "data/fastq_g_trimmed/{sample}.out",
        err=        "data/fastq_g_trimmed/{sample}.err",
        trimlog=    "data/fastq_g_trimmed/{sample}.log.gz"
    shell:
        """
        mkdir -p {params.dir}
        
        {trimmomatic} PE                            \
            -threads {threads}                      \
            -phred33                                \
            -trimlog >( gzip -9 > {log.trimlog} )   \
            {input.f}                               \
            {input.r}                               \
            {output.f}                              \
            {params.u1}                             \
            {output.r}                              \
            {params.u2}                             \
            HEADCROP:13                             \
            ILLUMINACLIP:{params.adaptors}:2:30:10  \
            MINLEN:31                               \
            AVGQUAL:10                              \
            MINLEN:31                               \
            TRAILING:19                             \
            MINLEN:31                               \
            TOPHRED33                               \
        >   {log.out}                               \
        2>  {log.err}
        
        gzip -dc {params.u1} {params.u2}            |
        gzip -9 > {output.u}
            
        rm {params.u1} {params.u2}
        """



rule make_bowtie2_index:
    input:
        assembly=   "data/assembly/assembly.fasta"
    output:
        name=       "data/index/assembly"        
    log:
        out=        "data/index/index.out",
        err=        "data/index/index.err"
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
        >   {log.out}           \
        2>  {log.err}

        touch {output}
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
        
        (( echo -e '@RG\tID:{params.sample}\tLB:truseq_{params.sample}\tPL:Illumina\tSM:{params.sample}' &&  \
        bowtie2                     \
            --threads   {threads}   \
            --phred33               \
            --quiet                 \
            --sensitive-local       \
            --no-unal               \
            -k  2                   \
            -x  {input.index}       \
            -1  {input.f}           \
            -2  {input.r}           \
            -U  {input.u}         ) |
        {samtools}  view            \
            -S  `# Input is SAM`    \
            -h  `# Print SAM header`\
            -u  `# Uncompressed SAM`\
            -                       |
        {samtools}  sort            \
            -@  {threads}           \
            -l  0                   \
            -o                      \
            - -                     |
        {samtools}  rmdup           \
            -s                      \
            - -                     |
        {samtools}  sort            \
            -@  {threads}           \
            -l  9                   \
            -o                      \
            - -                   ) \
        >   {output.bam}            \
        2>  {log.err}

        touch {log.out}
        """
    
