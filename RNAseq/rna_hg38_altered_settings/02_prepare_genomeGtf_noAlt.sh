#genome_noAlt.gtf is produced in order to run featureCounts generating results compatible with STAR counts output (primary for the purpose of comparing the results in order evaluate their similarity).

cd ./resources
sed '/_alt/d' genome.gtf > genome_noAlt.gtf
sed -i '/_decoy/d' genome_noAlt.gtf 
sed -i '/_random/d' genome_noAlt.gtf 
sed -i '/^HLA/d' genome_noAlt.gtf 
sed -i '/^chrUn_/d' genome_noAlt.gtf 

