#!/bin/bash
#
#SBATCH --job-name=snake_scRNA_seq
#
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=quake,normal
# send mail to this address
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuanxue@stanford.edu
#SBATCH --output=/oak/stanford/groups/quake/yuanxue/resources/sc_pipeline/snakemake_pipeline/singleCell_snake/log/snakemake.log.out
#SBATCH --error=/oak/stanford/groups/quake/yuanxue/resources/sc_pipeline/snakemake_pipeline/singleCell_snake/log/snakemake.log.err

config=$(basename $1)

MY_HOME=/oak/stanford/groups/quake/yuanxue
SCRIPTDIR=$MY_HOME/resources/sc_pipeline/snakemake_pipeline/singleCell_snake
SNAKEFILE=$SCRIPTDIR/snakefile
CONFIGFILE=$SCRIPTDIR/config/$1

# Specify log file
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE=$SCRIPTDIR/log/Snakefile.$DATETIME.log

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

echo Log is
echo $LOGFILE
echo $HOSTNAME
echo
