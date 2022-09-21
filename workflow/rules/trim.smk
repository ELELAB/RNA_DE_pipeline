rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    conda:
        "../wrappers/executive_wrappers/sra-tools/fasterq-dump/environment.yaml"
    script:
        "../wrappers/executive_wrappers/sra-tools/fasterq-dump/wrapper.py"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        "../../results/pipe/cutadapt/{sample}/{unit}.{fq}.{ext}",
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 1
    shell:
        "ln -s {input} {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
#        ["../../results/pipe/cutadapt/{sample}/{unit}.fq1.fastq.gz", "../../results/pipe/cutadapt/{sample}/{unit}.fq2.fastq.gz"]
    output:
        fastq1="../../results/trimmed/{sample}_{unit}_R1.fastq.gz",
        fastq2="../../results/trimmed/{sample}_{unit}_R2.fastq.gz",
        qc="../../results/trimmed/{sample}-{unit}_paired.qc.txt",
    log:
        "logs/cutadapt/{sample}_{unit}.log",
    params:
        extra=config["params"]["cutadapt-pe"],
        adapters="",
    threads: 19
    conda:
        "../wrappers/executive_wrappers/cutadapt/pe/environment.yaml"
    script:
        "../wrappers/executive_wrappers/cutadapt/pe/wrapper.py"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="../../results/trimmed/{sample}-{unit}_single.fastq.gz",
        qc="../../results/trimmed/{sample}-{unit}_single.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        extra=config["params"]["cutadapt-se"],
        adapters_r1=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 19
    conda:
        "../wrappers/executive_wrappers/cutadapt/se/environment.yaml"
    script:
        "../wrappers/executive_wrappers/cutadapt/se/wrapper.py"


rule merge_fastqs:
    input:
        get_fastqs,
    output:
        "../../results/merged/{sample}_{read}.fastq.gz",
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    shell:
        "cat {input} > {output} 2> {log}"
