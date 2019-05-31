#!/bin/bash
#
#SBATCH --job-name=snake_scRNA_seq
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=quake,normal
# send mail to this address
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuanxue@stanford.edu

config=$(basename $1)

MY_HOME=/oak/stanford/groups/quake/yuanxue
SCRIPTDIR=$MY_HOME/resources/sc_pipeline/snakemake_pipeline/singleCell_snake
SNAKEFILE=$SCRIPTDIR/snakefile
CONFIGFILE=$SCRIPTDIR/config/$1

# Specify parameters for snakemake
RESTART=1

if [ "$2" = "unlock" ]; then
    snakemake all -s $SNAKEFILE --configfile $CONFIGFILE --unlock
elif [ "$2" = 'dry' ]; then
    snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o log/{params.name}.%j.log" --keep-target-files -j 200 -w 120 -k --rerun-incomplete -n --restart-times $RESTART --quiet
elif [ "$2" = 'rerun' ]; then
    snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o log/{params.name}.%j.log" --keep-target-files -j 200 -w 120 -k --rerun-incomplete --restart-times $RESTART -F
else
    snakemake all --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o log/{params.name}.%j.log" --keep-target-files -j 200 -w 120 -k --rerun-incomplete --restart-times $RESTART
fi
