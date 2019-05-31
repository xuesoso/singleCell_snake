![Logo][scSnake.png]
# singleCell_snake
Snakemake pipeline for STAR alignment + htseq-count expression quantification for Smart-seq2 workflow.

Work in progress
----------------
+ Now: Add support for transcript to gene matrix transformation
+ Now: Add support for writing to anndata object
+ Now: Modularize rules
+ Now: Add sample dataset
+ Future: Enable support for chunky STAR alignment and htseq-count to save memory loading time.
+ Future: single-cell variant counting analysis

Prerequisites
-------------
+ STAR v2.60+
+ htseq-count v0.10.0+

Usage
-----
+ Modify parameters of yaml under config/

Replace {CUSTOM.YAML} with matching config.yaml under config/
Submit remote master job on cluster 
```bash
sbatch do.sh {CUSTOM.YAML}
```

To do a dry run
```bash
sbatch do.sh {CUSTOM.YAML} dry
```

To unlock previously failed snakemake run
```bash
sbatch do.sh {CUSTOM.YAML} unlock
```

To force a rerun
```bash
sbatch do.sh {CUSTOM.YAML} rerun
```
