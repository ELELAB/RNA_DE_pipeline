rule select_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/genotyped/all.vcf.gz",
    output:
        vcf=("results/filtered/all.{vartype}.vcf.gz"),
    params:
        extra=get_vartype_arg,
    log:
        "logs/gatk/selectvariants/{vartype}.log",
    wrapper:
        "0.74.0/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf="results/filtered/all.{vartype}.hardfiltered.vcf.gz",
    params:
        filters=get_filter,
    log:
        "logs/gatk/variantfiltration/{vartype}.log",
    wrapper:
        "0.74.0/bio/gatk/variantfiltration"


rule recalibrate_calls_snvs:
    input:
        vcf="results/filtered/all.snvs.vcf.gz",
        ref="resources/genome.fasta",
        hapmap="resources/edited_chr/hapmap_3.3.hg38.vcf.gz",
        omni="resources/edited_chr/1000G_omni2.5.hg38.vcf.gz",
        g1k="resources/edited_chr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="resources/edited_chr/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    output:
        vcf="results/filtered/all.snvs.recalibrated.vcf.gz",
        tranches="results/filtered/tranches/all.snvs.tranches"
    log:
        "logs/gatk/variantrecalibrator.snvs.log"
    params:
        mode="SNP",  # set mode, must be either SNP, INDEL or BOTH
        resources={"hapmap":{"known": False, "training": True, "truth": True, "prior": 15.0},
            "omni":{"known": False, "training": True, "truth": False, "prior": 12.0},
            "g1k":{"known": False, "training": True, "truth": False, "prior": 10.0},
            "dbsnp":{"known": True, "training": False, "truth": False, "prior": 2.0}},
#        annotation=[""],  # which fields to use with -an (see VariantRecalibrator docs)
#        annotation=["MappingQualityRankSumTest -an ReadPosRankSumTest -an FisherStrand -an MappingQuality -an StrandOddsRatio -an DepthPerAlleleBySample"],  # which fields to use with -an (see VariantRecalibrator docs)
        annotation=["QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP"],  #correct annotation format is based on input vcf (not well define in the literature)
        extra="--max-gaussians 6 --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0",  # optional
    wrapper:
        "0.74.0/bio/gatk/variantrecalibrator"

rule recalibrate_calls_indels:
    input:
        vcf="results/filtered/all.indels.vcf.gz",
        ref="resources/genome.fasta",
        mills="resources/edited_chr/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        axiomPoly="resources/edited_chr/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        dbsnp="resources/edited_chr/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
    output:
        vcf="results/filtered/all.indels.recalibrated.vcf.gz",
        tranches="results/filtered/tranches/all.indels.tranches"
    log:
        "logs/gatk/variantrecalibrator.indels.log"
    params:
        mode="INDEL",  # set mode, must be either SNP, INDEL or BOTH
        resources={"mills":{"known": False, "training": True, "truth": True, "prior": 12.0},
            "axiomPoly":{"known": False, "training": True, "truth": False, "prior": 10.0},
            "dbsnp":{"known": True, "training": False, "truth": False, "prior": 2.0}},
        annotation=["QD -an FS -an SOR -an DP"],
        extra="--max-gaussians 2 --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0",  # optional; originally 4
    wrapper:
        "0.74.0/bio/gatk/variantrecalibrator"

rule apply_vqsr_snvs:
    input:
        vcf="results/filtered/all.snvs.vcf.gz",
        recal="results/filtered/all.snvs.recalibrated.vcf.gz",
        tranches="results/filtered/tranches/all.snvs.tranches",
        ref="resources/genome.fasta",
    output:
        vcf="results/filtered/all.snvs_recal.vcf.gz"
    log:
        "logs/gatk/all.snvs.applyvqsr.log"
    params:
        mode="SNP",# set mode, must be either SNP, INDEL or BOTH
        extra="--truth-sensitivity-filter-level 99.0"  # optional
    resources:
        mem_mb=5000
    wrapper:
        "0.74.0/bio/gatk/applyvqsr"

rule apply_vqsr_indels:
    input:
        vcf="results/filtered/all.indels.vcf.gz",
        recal="results/filtered/all.indels.recalibrated.vcf.gz",
        tranches="results/filtered/tranches/all.indels.tranches",
        ref="resources/genome.fasta",
    output:
        vcf="results/filtered/all.indels_recal.vcf.gz"
    log:
        "logs/gatk/all.indels.applyvqsr.log"
    params:
        mode="INDEL",# set mode, must be either SNP, INDEL or BOTH
        extra="--truth-sensitivity-filter-level 99.0"  # optional
    resources:
        mem_mb=5000
    wrapper:
        "0.74.0/bio/gatk/applyvqsr"

rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{vartype}_{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype="recal"
            if config["filtering"]["vqsr"]
            else "hardfiltered",
        ),
    output:
        vcf="results/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
    wrapper:
        "0.74.0/bio/picard/mergevcfs"
