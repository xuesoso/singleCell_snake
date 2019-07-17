![Logo](scSnake.png)
# singleCell_snake
Snakemake pipeline for STAR alignment + htseq-count expression quantification for Smart-seq2 workflow.

Work in progress
----------------
~~+ Add support for transcript to gene matrix transformation~~

~~+ Modularize rules~~

~~+ Add sample dataset~~

~~+ single-cell variant counting analysis~~

+ Now: Add support for writing to anndata object
+ Now: Enable support for chunky STAR alignment and htseq-count to save memory loading time.

Prerequisites
-------------
+ STAR v2.60+
+ htseq-count v0.10.0+
+ Python 3+, Pandas, numpy
+ samtools v1.8+
+ bcftools v1.8+

Usage
-----
+ Modify parameters of yaml under config/

Replace {CUSTOM.YAML} with matching config.yaml under config/
Submit remote master job on cluster 
```bash
sh do.sh {CUSTOM.YAML}
```

To do a dry run
```bash
sh do.sh {CUSTOM.YAML} dry
```

To unlock previously failed snakemake run
```bash
sh do.sh {CUSTOM.YAML} unlock
```

To forceall on snakemake
```bash
sh do.sh {CUSTOM.YAML} forceall
```

To forcerun on snakemake
```bash
sh do.sh {CUSTOM.YAML} forcerun
```
