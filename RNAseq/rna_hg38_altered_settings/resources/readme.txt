Required resources:
1) genome.fasta
2) genome.gtf
For RNA seq, primary assembly should be used (without alt contigs, HLA, etc.)
In our case we create symbolic links to those files.

note: Resources were copied from Adrian's repo.

genome_noAlt.gtf is generated in order to run featureCounts producing results compatible with STAR counts output (primary for the purpose of comparing the results in order to evaluate pipelines' similarity).

