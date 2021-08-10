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
        "0.59.0/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf=("results/filtered/all.{vartype}.hardfiltered.vcf.gz"),
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
        hapmap="resources/nik_hapmap_3.3.hg38.vcf.gz",
        omni="resources/nik_1000G_omni2.5.hg38.vcf.gz",
        g1k="resources/nik_1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="resources/nik_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
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
#        annotation=["QualByDepth -an MappingQualityRankSumTest -an ReadPosRankSumTest -an FisherStrand -an RMSMappingQuality -an StrandOddsRatio -an Coverage"],
        annotation=["QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP"],
        extra="-AS --max-gaussians 6",  # optional
    wrapper:
        "0.74.0/bio/gatk/variantrecalibrator"

rule recalibrate_calls_indels:
    input:
        vcf="results/filtered/all.indels.vcf.gz",
        ref="resources/genome.fasta",
        hapmap="resources/nik_hapmap_3.3.hg38.vcf.gz",
        omni="resources/nik_1000G_omni2.5.hg38.vcf.gz",
        g1k="resources/nik_1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp="resources/nik_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    output:
        vcf="results/filtered/all.indels.recalibrated.vcf.gz",
	tranches="results/filtered/tranches/all.indels.tranches"
    log:
	"logs/gatk/variantrecalibrator.indels.log"
    params:
        mode="INDEL",  # set mode, must be either SNP, INDEL or BOTH
        resources={"hapmap":{"known": False, "training": True, "truth": True, "prior": 15.0},
	    "omni":{"known": False, "training": True, "truth": False, "prior": 12.0},
            "g1k":{"known": False, "training": True, "truth": False, "prior": 10.0},
            "dbsnp":{"known": True, "training": False, "truth": False, "prior": 2.0}},
#        annotation=[""],  # which fields to use with -an (see VariantRecalibrator docs)
#	annotation=["QualByDepth -an MappingQualityRankSumTest -an ReadPosRankSumTest -an FisherStrand -an StrandOddsRatio -an Coverage"],
        annotation=["QD -an FS -an SOR -an DP"],
        extra="-AS --max-gaussians 2",  # optional
    wrapper:
        "0.74.0/bio/gatk/variantrecalibrator"

rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{vartype}.{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype="recalibrated"
            if config["filtering"]["vqsr"]
            else "hardfiltered",
        ),
    output:
        vcf="results/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
    wrapper:
        "0.74.0/bio/picard/mergevcfs"
