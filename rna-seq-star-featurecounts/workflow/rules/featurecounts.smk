##################PROBLEM WITH DUPLICATED FILES BECAUSE OF 2 WILDCARDS???
rule aggregate_bam:
    input:
        lambda wc: get_star_output_all_units(wc, fi="bam"),
           expand(
               "results/star/pe/{unit.sample_name}-{unit.unit_name}/Aligned.out.bam",
                unit=units.itertuples(),
                ),
    output:
        "results/featurecounts/aggregated.txt"
    shell:
        "echo {input} > {output}"

###the same results as previous rule
#rule aggregate_bam:
#    input:
#        expand("results/star/pe/{sample}-{unit}/Aligned.out.bam", sample=units.sample_name, unit=units.unit_name)
#    output:
#        "results/featurecounts/aggregated.txt"
#    shell:
#        "echo {input} > {output}"

###"lane1" is hardcoded for now (it's necessary to remove duplicates from the list created by "expand" function
rule feature_counts:
    input:
        sam=expand("results/star/pe/{sample}-lane1/Aligned.out.bam", sample=units.sample_name, unit=units.unit_name),
#        with open('results/featurecounts/aggregated.txt','r') as file:
#            sam = file.read()
#            print(sam)
#print "results/featurecounts/aggregated.txt"
        annotation="resources/genome.gtf",
        # optional input
        # chr_names="",          # implicitly sets the -A flag
        fasta="resources/genome.fasta" # implicitly sets the -G flag,

    output:
        "results/featurecounts/all.featureCounts",
        "results/featurecounts/all.featureCounts.summary",
        "results/featurecounts/all.featureCounts.jcounts", 

    threads:
        10
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2"
    log:
        "logs/featurecounts/B.log"
    wrapper:
        "v0.75.0/bio/subread/featurecounts"

