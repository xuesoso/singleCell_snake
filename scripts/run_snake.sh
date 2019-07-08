#!/bin/bash
#
#SBATCH --job-name=snake_scRNA_seq
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=owners,normal
# send mail to this address
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuanxue@stanford.edu

config=$(basename $1)

MY_HOME=/oak/stanford/groups/quake/yuanxue
SCRIPTDIR=$MY_HOME/resources/sc_pipeline/snakemake_pipeline/singleCell_snake
SNAKEFILE=$SCRIPTDIR/snakefile
SNAKEFILE_SNP=$SCRIPTDIR/snakefile_snp
CONFIGFILE=$SCRIPTDIR/config/$1
NJOBS=200
WAIT=120

# Specify parameters for snakemake
RESTART=1

# Load bcftools if on Stanford Sherlock, otherwise you need to make sure system has access to bcftools
module load biology
module load bcftools/1.8
#

if [ "$2" = "unlock" ]; then
    snakemake all -s $SNAKEFILE --configfile $CONFIGFILE --unlock
elif [ "$2" = 'dry' ]; then
    snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete -n --restart-times $RESTART --quiet
elif [ "$2" = 'forceall' ]; then
    snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART -F
elif [ "$2" = 'forcerun' ]; then
    snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART -R
elif [ "$2" = 'snp' ]; then
    snakemake all --snakefile $SNAKEFILE_SNP --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART
else
    snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART
fi
