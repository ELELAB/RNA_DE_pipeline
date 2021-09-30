#rule get_genome:
#    output:
#        "resources/genome.fasta",
#    log:
#        "logs/get-genome.log",
#    params:
#        species=config["ref"]["species"],
#        datatype="dna",
#        build=config["ref"]["build"],
#        release=config["ref"]["release"],
#    cache: True
#    conda:
#        "../wrappers/executive_wrappers/reference/ensembl-sequence/environment.yaml"
#    script:
#        "../wrappers/executive_wrappers/reference/ensembl-sequence/wrapper.py"

#rule get_annotation:
#    output:
#        "resources/genome.gtf",
#    params:
#        species=config["ref"]["species"],
#        fmt="gtf",
#        build=config["ref"]["build"],
#        release=config["ref"]["release"],
#        flavor="",
#    cache: True
#    log:
#        "logs/get_annotation.log",
#    conda:
#        "../wrappers/executive_wrappers/reference/ensembl-annotation/environment.yaml"
#    script:
#        "../wrappers/executive_wrappers/reference/ensembl-annotation/wrapper.py"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    conda:
        "../wrappers/executive_wrappers/samtools/index/environment.yaml"
    script:
        "../wrappers/executive_wrappers/samtools/index/wrapper.py"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    conda:
        "../wrappers/executive_wrappers/bwa/index/environment.yaml"
    script:
        "../wrappers/executive_wrappers/bwa/index/wrapper.py"


rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 4
    params:
        extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 149",
    log:
        "logs/star_index_genome.log",
    cache: True
    conda:
        "../wrappers/executive_wrappers/star/index/environment.yaml"
    script:
        "../wrappers/executive_wrappers/star/index/wrapper.py"
