rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        index="resources/star_genome",
    output:
        "../../results/star/pe/{sample}-{unit}/{sample}_Aligned.out.bam",
        "../../results/star/pe/{sample}-{unit}/{sample}_ReadsPerGene.out.tab",
    log:
        "logs/star-pe/{sample}-{unit}.log",
    params:
        FileNamePrefix="../../results/star/pe/{sample}-{unit}/{sample}_",
        index=lambda wc, input: input.index,
        extra="--sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"],
        ),
    threads: 6
    conda:
        "../wrappers/executive_wrappers/star/align/environment.yaml"
    script:
        "../wrappers/executive_wrappers/star/align/wrapper_nik.py"


rule picard_sort_sam_pe:
    input:
        bam="../../results/star/pe/{sample}-{unit}/{sample}_Aligned.out.bam",
    output:
        "../../results/star/pe/{sample}-{unit}/{sample}_Aligned.out.sorted.bam"
    params:
        sort_order="coordinate", 
        extra="VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
    log:
        "logs/picard/sortSam/{sample}-{unit}.log"
    resources:
        mem_mb=10024
    conda:
        "../wrappers/executive_wrappers/picard/sortsam/environment.yaml"
    script:
        "../wrappers/executive_wrappers/picard/sortsam/wrapper.py"


rule align_se:
    input:
        fq1=get_map_reads_input_R1,
        index="resources/star_genome",
    output:
        "../../results/star/se/{sample}-{unit}/{sample}_Aligned.out.bam",
        "../../results/star/se/{sample}-{unit}/{sample}_ReadsPerGene.out.tab",
    log:
        "logs/star-se/{sample}-{unit}.log",
    params:
        FileNamePrefix="../../results/star/pe/{sample}-{unit}/{sample}_",
        index=lambda wc, input: input.index,
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: 19
    conda:
        "../wrappers/executive_wrappers/star/align/environment.yaml"
    script:
        "../wrappers/executive_wrappers/star/align/wrapper.py"


rule picard_sort_sam_se:
    input:
        bam="../../results/star/se/{sample}-{unit}/{sample}_Aligned.out.bam",
    output:
        "../../results/star/se/{sample}-{unit}/{sample}_Aligned.out.sorted.bam"
    params:
        sort_order="coordinate", 
        extra="VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"
    log:
        "logs/picard/sortSam/{sample}-{unit}.log"
    resources:
        mem_mb=10024
    conda:
        "../wrappers/executive_wrappers/picard/sortsam/environment.yaml"
    script:
        "../wrappers/executive_wrappers/picard/sortsam/wrapper.py"
