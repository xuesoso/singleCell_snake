# Specify log file
MY_HOME=/oak/stanford/groups/quake/yuanxue
SCRIPTDIR=$MY_HOME/resources/sc_pipeline/snakemake_pipeline/singleCell_snake
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE_OUT=$SCRIPTDIR/log/Snakefile.$DATETIME.$1.log.out
LOGFILE_ERR=$SCRIPTDIR/log/Snakefile.$DATETIME.$1.log.err

sbatch --output=$LOGFILE_OUT --error=$LOGFILE_ERR $SCRIPTDIR/scripts/run_snake.sh $1 $2

echo Log is
echo $LOGFILE_OUT
echo $HOSTNAME
