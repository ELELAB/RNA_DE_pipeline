## FASTQC and RSEQC
rule fastqc_input:
    input:
        prepare_fastqc_input,
    output:
        "../../../../results/fastqc_input/{sample}-{unit}-{fq}.{ext}",
    log:
        "logs/fastqc_input/{sample}-{unit}-{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0
    shell:
        "ln -s {input} {output} 2> {log}"


rule fastqc_qc:
    input:
        "../../../../results/fastqc_input/{sample}-{unit}-{fq}.fastq.gz",
    output:
        html="../../../../qc/fastqc/{sample}-{unit}-{fq}_fastqc.html",
        zip="../../../../qc/fastqc/{sample}-{unit}-{fq}_fastqc.zip",
    priority: 1
    log:
        "logs/fastqc/{sample}-{unit}-{fq}.log",
    params:
        "",
    conda:
        "../wrappers/executive_wrappers/fastqc/environment.yaml"
    script:
        "../wrappers/executive_wrappers/fastqc/wrapper.py"


rule rseqc_gtf2bed:
    input:
        "resources/genome.gtf",
    output:
        bed="../../../../qc/rseqc/annotation.bed",
        db=temp("../../../../qc/rseqc/annotation.db"),
    log:
        "logs/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam=get_star_bam,
        bed="../../../../qc/rseqc/annotation.bed",
    output:
        "../../../../qc/rseqc/{sample}-{unit}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam=get_star_bam,
        bed="../../../../qc/rseqc/annotation.bed",
    output:
        "../../../../qc/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        get_star_bam,
    output:
        "../../../../qc/rseqc/{sample}-{unit}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam=get_star_bam,
        bed="../../../../qc/rseqc/annotation.bed",
    output:
        "../../../../qc/rseqc/{sample}-{unit}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam=get_star_bam,
        bed="../../../../qc/rseqc/annotation.bed",
    output:
        "../../../../qc/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam=get_star_bam,
        bed="../../../../qc/rseqc/annotation.bed",
    output:
        "../../../../qc/rseqc/{sample}-{unit}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        get_star_bam,
    output:
        "../../../../qc/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        get_star_bam,
    output:
        "../../../../qc/rseqc/{sample}-{unit}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}-{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


#In case of Cutadapt (also other tools), Multiqc is parsing file names from *.qc.txt. This cause the problem because cutadapt input is saved as "sample per folder". That results in overwriting an input from Cutadapt by Multiqc and reporting only one sample, eg. "lane1". 
rule multiqc:
    input:
        lambda wc: get_star_output_all_units(wc, fi="bam"),
        expand(
            "../../../../qc/fastqc/{unit.sample_name}-{unit.unit_name}-fq1_fastqc.zip",
            unit=units.itertuples(),  
        ),
        expand(
            "../../../../qc/fastqc/{unit.sample_name}-{unit.unit_name}-fq2_fastqc.zip",
            unit=units.itertuples(),  
        ),
#        expand(
#            "../../../../results/trimmed/{unit.sample_name}-{unit.unit_name}_paired.qc.txt",
#            unit=units.itertuples(),
#	), 
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.junctionanno.junction.bed",
            unit=units.itertuples(),
        ),
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.junctionsat.junctionSaturation_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.infer_experiment.txt",
            unit=units.itertuples(),
        ),
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.stats.txt",
            unit=units.itertuples(),
        ),
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.inner_distance_freq.inner_distance.txt",
            unit=units.itertuples(),
        ),
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.readdistribution.txt",
            unit=units.itertuples(),
        ),
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.readdup.DupRate_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "../../../../qc/rseqc/{unit.sample_name}-{unit.unit_name}.readgc.GC_plot.pdf",
            unit=units.itertuples(),
        ),
        expand(
            "logs/rseqc/rseqc_junction_annotation/{unit.sample_name}-{unit.unit_name}.log",
            unit=units.itertuples(),
        ),
    output:
        "../../../../qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../wrappers/executive_wrappers/multiqc/environment.yaml"
    script:
        "../wrappers/executive_wrappers/multiqc/wrapper.py"

