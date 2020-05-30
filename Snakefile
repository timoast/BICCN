"""
Process scATAC-seq data
Requires bwa-mem, samtools, sinto, bgzip, and tabix
"""

IND, = glob_wildcards("regions/{rep}.txt")


rule all:
    input: expand("fragments/{rep}.sort.bed.gz", rep=IND)

rule get_genome:
    """Download genome and build bwa index"""
    output:
        "genome/mm10.fa.gz"
    threads: 1
    shell:
        """
        cd genome
        wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
        """

rule bwa_build:
    """Build index for bwa mem"""
    input:
        "genome/mm10.fa.gz"
    output:
        "genome/mm10.fa"
    threads: 1
    shell:
        """
        gzip -d genome/mm10.fa.gz
        bwa index genome/mm10.fa
        """

rule download:
    """
    Download fastq files for each replicate
    expand tar
    cat fastqs from different runs
    """
    input:
        "regions/{rep}.txt"
    output:
        "fastq/{rep}/read1.fastq.gz"
    threads: 1
    shell:
        """
        wget -i {input} -P fastq/{wildcards.rep}
        cd fastq/{wildcards.rep}
        for i in *.tar; do
            tar -xf $i
        done
        mv */* .
        cat *R1.fastq.gz > read1.fastq.gz
        cat *R2.fastq.gz > read2.fastq.gz
        rm *R1.fastq.gz *R2.fastq.gz
        """

rule align:
    input:
        genome = "genome/mm10.fa",
        reads = "fastq/{rep}/read1.fastq.gz"
    output:
        "mapped/{rep}.sort.bam"
    threads: 10
    shell:
        """
        bwa mem -t {threads} {input.genome} \
            fastq/{wildcards.rep}/read1.fastq.gz \
            fastq/{wildcards.rep}/read2.fastq.gz \
            | samtools view -b - > fastq/{wildcards.rep}/aln.bam
        
        samtools sort -@ {threads} fastq/{wildcards.rep}/aln.bam -o {output}
        samtools index -@ {threads} {output}

        rm fastq/{wildcards.rep}/aln.bam
        """

rule fragments:
    input:
        "mapped/{rep}.sort.bam"
    output:
        "fragments/{rep}.bed"
    threads: 10
    shell:
        """
        sinto fragments -b {input} -p {threads} -f fragments/{wildcards.rep}.bed --barcode_regex "[^:]*"
        """

rule sort_frags:
    input:
        "fragments/{rep}.bed"
    output:
        "fragments/{rep}.sort.bed.gz"
    threads: 5
    shell:
        """
        sort -k1,1 -k2,2n fragments/{wildcards.rep}.bed > fragments/{wildcards.rep}.sort.bed
        bgzip -@ {threads} fragments/{wildcards.rep}.sort.bed
        tabix -p bed fragments/{wildcards.rep}.sort.bed.gz
        rm fragments/{wildcards.rep}.bed
        """