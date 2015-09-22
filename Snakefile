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
samtools=       "./bin/samtools" # v0.1.20
bcftools=       "./bin/bcftools" # v0.1.20

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
            sample= gsamples),
        "data/g2t_k2/all.bam",
        "data/g2t_k2/all.bam.bai",
        "data/g2t_k2/all.bcf",
        "data/g2t_k2/all.vcf",
        #expand("data/g2t_k1/{sample}.bam",
        #    sample= gsamples),
        #"data/g2t_k1/all.bam",
        #"data/g2t_k1/all.bcf",
        #"data/g2t_k1/all.vcf"



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



rule index_assembly:
    input:
        assembly=   "data/assembly/assembly.fasta"
    output:
        index=      "data/assembly/assembly.fasta.fai"
    log:
        out=        "data/assembly/assembly_faidx.out",
        err=        "data/assembly/assembly_faidx.err"
    threads:
        1
    shell:
        """
        {samtools} faidx {input.assembly}   \
        >   {log.out}                       \
        2>  {log.err}
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



#############
## G2T k=2 ##
#############

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



rule g2t_k2_merge:
    input:
        bams=       expand("data/g2t_k2/{samples}.bam",
                        samples= gsamples)
    output:
        bam=        "data/g2t_k2/all.bam"
    log:
        out=        "data/g2t_k2/all.out",
        err=        "data/g2t_k2/all.err"
    params:
        bam_tmp=    "data/g2t_k2/all_tmp.bam",
        header=     "data/g2t_k2/header.txt",
        header1=    "data/g2t_k2/tmp1",
        header2=    "data/g2t_k2/tmp2",
        header3=    "data/g2t_k2/tmp3"
    threads:
        24
    shell:
        """
        {samtools} merge        \
            -f                  \
            -r                  \
            -l  9               \
            -@  {threads}       \
            {params.bam_tmp}    \
            {input.bams}        \
        >   {log.out}           \
        2>  {log.err}

        {samtools} view         \
            -H                  \
            {params.bam_tmp}    \
        >   {params.header}     \
        2>> {log.err}
        
        head    -n  -1  {params.header} >   {params.header1}
        tail    -n  1   {params.header} >   {params.header3}
        
        cat /dev/null > {params.header2}
        
        for sample in {gsamples}; do
            echo -e "@RG\tID:${{sample}}\tLB:truseq_${{sample}}\tPL:Illumina\tSM:${{sample}}" >> {params.header2}
        done
        
        cat {params.header1} {params.header2} {params.header3} > {params.header}

        {samtools} reheader     \
            {params.header}     \
            {params.bam_tmp}    \
        >   {output.bam}        \
        2>> {log.err}

        rm  {params}
        """



rule g2t_k2_bai:
    input:
        bam=    "data/g2t_k2/all.bam"
    output:
        bai=    "data/g2t_k2/all.bam.bai"
    log:
        out=    "data/g2t_k2/all_index.out",
        err=    "data/g2t_k2/all_index.err"
    threads:
        1
    shell:
        """
        {samtools}  index   \
            {input.bam}     \
        >   {log.out}       \
        2>  {log.err}
        """



rule g2t_k2_bcf:
    """
    From the BAM and the assembly, call snps with samtools mpileup and generate 
    a BCF file.
    Explanation of all the flags:
        Input options:
            -A: Count anomalous flags
            -B: Disable BAQ computation
            -f: faidx indexed reference sequence file [null]
            -C: parameter for adjusting mapQ; 0 to disable [0]
        Output options:
            -D: output per-sample DP in BCF (require -g/-u)
            -g: generate BCF output (genotype likelihoods)
            -u: generate uncompress BCF output
        SNP/INDEL genotype likelyhoods options
            -I: do not perform indel calling
    """
    input:
        assembly=       "data/assembly/assembly.fasta",
        assembly_idx=   "data/assembly/assembly.fasta.fai",
        bam=            "data/g2t_k2/all.bam"
    output:
        bcf=            "data/g2t_k2/all.bcf"
    threads:
        1
    log:
        out=            "data/g2t_k2/all_call_snps.out",
        err=            "data/g2t_k2/all_call_snps.err"
    params:
        mapq= "50"
    shell:
        """
        {samtools}  mpileup         \
            -A                      \
            -B                      \
            -D                      \
            -g                      \
            `#-u`                   \
            -I                      \
            -f  {input.assembly}    \
            -C  {params.mapq}       \
            {input.bam}             \
        >   {output.bcf}            \
        2>  {log.err}

        touch {log.out}
        """



rule g2t_k2_create_samples:
    output:
        samples=    "data/g2t_k2/samples.txt"
    threads:
        1
    log:
        out=        "data/g2t_k2/samples.out",
        err=        "data/g2t_k2/samples.err"
    shell:
        """
        cat /dev/null > {output.samples}

        for sample in {gsamples}; do
            echo -e "${{sample}}\t2" >> {output.samples}
        done
        """



rule g2t_k2_vcf:
    """
        Input/output options:
            -L: calculate LD for adjacent sites
            -N: skip sites where REF is not A/C/G/T
        Consensus/variant calling options:
            -c: SNP calling (force -e)
            -e: likelihood based analyses
            -g: call genotypes at variant sites (force -c)
            -I: skip indels
            -v: output potential variant sites only (force -c)
    """
    input:
        samples=    "data/g2t_k2/samples.txt",
        bcf=        "data/g2t_k2/all.bcf"
    output:
        vcf=        "data/g2t_k2/all.vcf"
    log:
        out=        "data/g2t_k2/vcf.out",
        err=        "data/g2t_k2/vcf.err"
    threads:
        1
    shell:
        """
        {bcftools} view         \
            -LNcegIv            \
            -s  {input.samples} \
            {input.bcf}         \
        >   {output.vcf}        \
        2>  {log.err}

        touch   {log.out} 
        """



#############
## G2T k=1 ##
#############



rule make_bowtie2_index:
    input:
        assembly=   "data/assembly/filtered.fasta"
    output:
        name=       "data/index/filtered"        
    log:
        out=        "data/index/filtered.out",
        err=        "data/index/filtered.err"
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



rule g2t_k1_mapping:
    """
    Perform the g2t mapping with bowtie2 and the 
    sam | bam | rmdup | sort conversions
    """
    input:
        index=  "data/index/filtered",
        f=      "data/fastq_g_trimmed/{sample}_1.fastq.gz",
        r=      "data/fastq_g_trimmed/{sample}_2.fastq.gz",
        u=      "data/fastq_g_trimmed/{sample}_U.fastq.gz"
    output:
        bam=    "data/g2t_k1/{sample}.bam"
    log:
        out=    "data/g2t_k1/{sample}.out",
        err=    "data/g2t_k1/{sample}.err"
    params:
        sample= "{sample}",
        dir=    "data/g2t_k1"
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
            -k  1                   \
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


rule g2t_k1_merge:
    input:
        bams=       expand("data/g2t_k1/{samples}.bam",
                        samples= gsamples)
    output:
        bam=        "data/g2t_k1/all.bam"
    log:
        out=        "data/g2t_k1/all.out",
        err=        "data/g2t_k1/all.err"
    params:
        bam_tmp=    "data/g2t_k1/all_tmp.bam",
        header=     "data/g2t_k1/header.txt",
        header1=    "data/g2t_k1/tmp1",
        header2=    "data/g2t_k1/tmp2",
        header3=    "data/g2t_k1/tmp3"
    threads:
        24
    shell:
        """
        {samtools} merge        \
            -f                  \
            -r                  \
            -l  9               \
            -@  {threads}       \
            {params.bam_tmp}    \
            {input.bams}        \
        >   {log.out}           \
        2>  {log.err}

        {samtools} view         \
            -H                  \
            {params.bam_tmp}    \
        >   {params.header}     \
        2>> {log.err}
        
        head    -n  -1  {params.header} >   {params.header1}
        tail    -n  1   {params.header} >   {params.header3}
        
        cat /dev/null > {params.header2}
        
        for sample in {gsamples}; do
            echo -e "@RG\tID:${{sample}}\tLB:truseq_${{sample}}\tPL:Illumina\tSM:${{sample}}" >> {params.header2}
        done
        
        cat {params.header1} {params.header2} {params.header3} > {params.header}

        {samtools} reheader     \
            {params.header}     \
            {params.bam_tmp}    \
        >   {output.bam}        \
        2>> {log.err}

        rm  {params}
        """



rule g2t_k1_bai:
    input:
        bam=    "data/g2t_k1/all.bam"
    output:
        bai=    "data/g2t_k1/all.bam.bai"
    log:
        out=    "data/g2t_k1/all_index.out",
        err=    "data/g2t_k1/all_index.err"
    threads:
        1
    shell:
        """
        {samtools}  index   \
            {input.bam}     \
        >   {log.out}       \
        2>  {log.err}
        """



rule g2t_k1_bcf:
    """
    From the BAM and the assembly, call snps with samtools mpileup and generate 
    a BCF file.
    Explanation of all the flags:
        Input options:
            -A: Count anomalous flags
            -B: Disable BAQ computation
            -f: faidx indexed reference sequence file [null]
            -C: parameter for adjusting mapQ; 0 to disable [0]
        Output options:
            -D: output per-sample DP in BCF (require -g/-u)
            -g: generate BCF output (genotype likelihoods)
            -u: generate uncompress BCF output
        SNP/INDEL genotype likelyhoods options
            -I: do not perform indel calling
    """
    input:
        assembly=       "data/assembly/assembly.fasta",
        assembly_idx=   "data/assembly/assembly.fasta.fai",
        bam=            "data/g2t_k1/all.bam"
    output:
        bcf=            "data/g2t_k1/all.bcf"
    threads:
        1
    log:
        out=            "data/g2t_k1/all_call_snps.out",
        err=            "data/g2t_k1/all_call_snps.err"
    params:
        mapq= "50"
    shell:
        """
        {samtools}  mpileup         \
            -A                      \
            -B                      \
            -D                      \
            -g                      \
            `#-u`                   \
            -I                      \
            -f  {input.assembly}    \
            -C  {params.mapq}       \
            {input.bam}             \
        >   {output.bcf}            \
        2>  {log.err}

        touch {log.out}
        """



rule g2t_k1_create_samples:
    output:
        samples=    "data/g2t_k1/samples.txt"
    threads:
        1
    log:
        out=        "data/g2t_k1/samples.out",
        err=        "data/g2t_k1/samples.err"
    shell:
        """
        cat /dev/null > {output.samples}

        for sample in {gsamples}; do
            echo -e "${{sample}}\t2" >> {output.samples}
        done
        """



rule g2t_k1_vcf:
    """
        Input/output options:
            -L: calculate LD for adjacent sites
            -N: skip sites where REF is not A/C/G/T
        Consensus/variant calling options:
            -c: SNP calling (force -e)
            -e: likelihood based analyses
            -g: call genotypes at variant sites (force -c)
            -I: skip indels
            -v: output potential variant sites only (force -c)
    """
    input:
        samples=    "data/g2t_k1/samples.txt",
        bcf=        "data/g2t_k1/all.bcf"
    output:
        vcf=        "data/g2t_k1/all.vcf"
    log:
        out=        "data/g2t_k1/vcf.out",
        err=        "data/g2t_k1/vcf.err"
    threads:
        1
    shell:
        """
        {bcftools} view         \
            -LNcegIv            \
            -s  {input.samples} \
            {input.bcf}         \
        >   {output.vcf}        \
        2>  {log.err}

        touch   {log.out} 
        """



################################################################################
## T2T                                                                        ##
################################################################################

rule t2t_mapping:
    """
    Perform the t2t mapping with bowtie2 and the 
    sam | bam | rmdup | sort conversions
    """
    input:
        index=  "data/index/filtered",
        f=      "data/fastq_t_trimmed/{sample}_1.fastq.gz",
        r=      "data/fastq_t_trimmed/{sample}_2.fastq.gz",
        u=      "data/fastq_t_trimmed/{sample}_U.fastq.gz"
    output:
        bam=    "data/t2t/{sample}.bam"
    log:
        out=    "data/t2t/{sample}.out",
        err=    "data/t2t/{sample}.err"
    params:
        sample= "{sample}",
        dir=    "data/t2t"
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



rule t2t_merge:
    input:
        bams=       expand("data/t2t/{samples}.bam",
                        samples= tsamples)
    output:
        bam=        "data/t2t/all.bam"
    log:
        out=        "data/t2t/all.out",
        err=        "data/t2t/all.err"
    params:
        bam_tmp=    "data/t2t/all_tmp.bam",
        header=     "data/t2t/header.txt",
        header1=    "data/t2t/tmp1",
        header2=    "data/t2t/tmp2",
        header3=    "data/t2t/tmp3"
    threads:
        24
    shell:
        """
        {samtools} merge        \
            -f                  \
            -r                  \
            -l  9               \
            -@  {threads}       \
            {params.bam_tmp}    \
            {input.bams}        \
        >   {log.out}           \
        2>  {log.err}

        {samtools} view         \
            -H                  \
            {params.bam_tmp}    \
        >   {params.header}     \
        2>> {log.err}
        
        head    -n  -1  {params.header} >   {params.header1}
        tail    -n  1   {params.header} >   {params.header3}
        
        cat /dev/null > {params.header2}
        
        for sample in {gsamples}; do
            echo -e "@RG\tID:${{sample}}\tLB:truseq_${{sample}}\tPL:Illumina\tSM:${{sample}}" >> {params.header2}
        done
        
        cat {params.header1} {params.header2} {params.header3} > {params.header}

        {samtools} reheader     \
            {params.header}     \
            {params.bam_tmp}    \
        >   {output.bam}        \
        2>> {log.err}

        rm  {params}
        """



rule t2t_bai:
    input:
        bam=    "data/t2t/all.bam"
    output:
        bai=    "data/t2t/all.bam.bai"
    log:
        out=    "data/t2t/all_index.out",
        err=    "data/t2t/all_index.err"
    threads:
        1
    shell:
        """
        {samtools}  index   \
            {input.bam}     \
        >   {log.out}       \
        2>  {log.err}
        """



rule t2t_bcf:
    """
    From the BAM and the assembly, call snps with samtools mpileup and generate 
    a BCF file.
    Explanation of all the flags:
        Input options:
            -A: Count anomalous flags
            -B: Disable BAQ computation
            -f: faidx indexed reference sequence file [null]
            -C: parameter for adjusting mapQ; 0 to disable [0]
        Output options:
            -D: output per-sample DP in BCF (require -g/-u)
            -g: generate BCF output (genotype likelihoods)
            -u: generate uncompress BCF output
        SNP/INDEL genotype likelyhoods options
            -I: do not perform indel calling
    """
    input:
        assembly=       "data/assembly/filtered.fasta",
        assembly_idx=   "data/assembly/filtered.fasta.fai",
        bam=            "data/t2t/all.bam"
    output:
        bcf=            "data/t2t/all.bcf"
    threads:
        1
    log:
        out=            "data/t2t/all_call_snps.out",
        err=            "data/t2t/all_call_snps.err"
    params:
        mapq= "50"
    shell:
        """
        {samtools}  mpileup         \
            -A                      \
            -B                      \
            -D                      \
            -g                      \
            `#-u`                   \
            -I                      \
            -f  {input.assembly}    \
            -C  {params.mapq}       \
            {input.bam}             \
        >   {output.bcf}            \
        2>  {log.err}

        touch {log.out}
        """



rule t2t_create_samples:
    output:
        samples=    "data/t2t/samples.txt"
    threads:
        1
    log:
        out=        "data/t2t/samples.out",
        err=        "data/t2t/samples.err"
    shell:
        """
        cat /dev/null > {output.samples}

        for sample in {gsamples}; do
            echo -e "${{sample}}\t2" >> {output.samples}
        done
        """



rule g2t_k2_vcf:
    """
        Input/output options:
            -L: calculate LD for adjacent sites
            -N: skip sites where REF is not A/C/G/T
        Consensus/variant calling options:
            -c: SNP calling (force -e)
            -e: likelihood based analyses
            -g: call genotypes at variant sites (force -c)
            -I: skip indels
            -v: output potential variant sites only (force -c)
    """
    input:
        samples=    "data/g2t_k2/samples.txt",
        bcf=        "data/g2t_k2/all.bcf"
    output:
        vcf=        "data/g2t_k2/all.vcf"
    log:
        out=        "data/g2t_k2/vcf.out",
        err=        "data/g2t_k2/vcf.err"
    threads:
        1
    shell:
        """
        {bcftools} view         \
            -LNcegIv            \
            -s  {input.samples} \
            {input.bcf}         \
        >   {output.vcf}        \
        2>  {log.err}

        touch   {log.out} 
        """
       
