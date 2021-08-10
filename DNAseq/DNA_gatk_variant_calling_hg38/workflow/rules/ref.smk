rule get_genome:
    output:
        "resources/vep/genome.fasta",
    log:
        "logs/vep/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    script:
        "../wrappers/ensembl-sequence_0.74.0/script.py"


checkpoint genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/vep/genome.fasta",
    output:
        "resources/vep/genome.dict",
    log:
        "logs/vep/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/vep/genome.fasta.fai",
    output:
        vcf="resources/vep/variation.vcf.gz",
    log:
        "logs/vep/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "resources/vep/variation.vcf.gz",
    output:
        "resources/vep/variation.noiupac.vcf.gz",
    log:
        "logs/vep/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: False
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants:
    input:
        "resources/vep/variation.noiupac.vcf.gz",
    output:
        "resources/vep/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/vep/tabix/variation.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "0.74.0/bio/tabix"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa", ".alt"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.74.0/bio/bwa/index"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "0.74.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "0.74.0/bio/vep/plugins"
