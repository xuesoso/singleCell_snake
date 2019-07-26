#!/bin/bash
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# send mail to this address
#SBATCH --mail-type=END,FAIL

config=$(basename $1)

MY_HOME=/oak/stanford/groups/quake/yuanxue
SCRIPTDIR=$MY_HOME/resources/sc_pipeline/snakemake_pipeline/singleCell_snake
SNAKEFILE=$SCRIPTDIR/snakefile
SNP_RULE=bam_to_vcf
CONFIGFILE=$SCRIPTDIR/config/$1
if [ "$2" = "local" ]; then
    NJOBS=10
else
    NJOBS=200
    # Load bcftools if on Stanford Sherlock, otherwise
    # you need to make sure system has access to bcftools & samtools
    module load biology
    module load bcftools/1.8
    module load samtools
fi

if [ -z "$3" ]; then
    echo "Running rule all"
    TARGETRULE=all
else    
    echo "Running rule $3"
    TARGETRULE=$3
fi
WAIT=120

# Specify parameters for snakemake
RESTART=1

if [ "$2" = "unlock" ]; then
    snakemake $TARGETRULE -s $SNAKEFILE --configfile $CONFIGFILE --unlock
elif [ "$2" = 'dry' ]; then
    snakemake $TARGETRULE --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete -n --restart-times $RESTART --quiet
elif [ "$2" = 'forceall' ]; then
    snakemake $TARGETRULE --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART -F
elif [ "$2" = 'forcerun' ]; then
    snakemake $TARGETRULE --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART -R
elif [ "$2" = 'local' ]; then
    snakemake $TARGETRULE --snakefile $SNAKEFILE --configfile $CONFIGFILE --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART
else
    snakemake $TARGETRULE --snakefile $SNAKEFILE --configfile $CONFIGFILE --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} --time={params.time} -o log/{params.name}.%j.log" --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART
fi
