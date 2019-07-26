#!/bin/bash

# Users: Please specify cluster partition if running on cluster 
SCRIPTDIR="$( cd "$(dirname "$0")" ; pwd -P )"
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
LOGFILE_OUT=$SCRIPTDIR/log/Snakefile.$DATETIME.$1.log.out
LOGFILE_ERR=$SCRIPTDIR/log/Snakefile.$DATETIME.$1.log.err
PARTITION=owners,normal,quake
JOBNAME=snake_scRNA_seq

if [ "$2" = "local" ]; then
    sh $SCRIPTDIR/scripts/run_snake.sh $1 $2 $3
else
    sbatch --job-name=$JOBNAME --partition=$PARTITION --output=$LOGFILE_OUT --error=$LOGFILE_ERR $SCRIPTDIR/scripts/run_snake.sh $1 $2 $3
fi

echo Log is
echo $LOGFILE_OUT
echo $HOSTNAME
