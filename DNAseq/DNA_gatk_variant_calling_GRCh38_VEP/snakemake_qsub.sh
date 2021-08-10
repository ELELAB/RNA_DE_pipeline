#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=dtu_00011 -A dtu_00011
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N nik_snakemake
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e test.err
#PBS -o test.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=15
### Memory
#PBS -l mem=250gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 12 hours)
#PBS -l walltime=60:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
MY_dir=/home/projects/dtu_00011/people/niktom/snakemake/test_gatk_local2
echo Working directory is $MY_dir
cd $MY_dir

### Here follows the user commands:
# Define number of processors
NPROCS=30
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load tools
module load miniconda3/4.10.1
conda init --all
conda activate /home/people/niktom/.conda/envs/snakemake/

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here
snakemake --cores 20 --use-conda --rerun-incomplete
