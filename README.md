![Logo](scSnake.png)
# singleCell_snake
Snakemake pipeline for STAR alignment + htseq-count expression quantification for Smart-seq2 workflow.

Work in progress
----------------
~~+ Now: Add support for transcript to gene matrix transformation~~

~~+ Now: Modularize rules~~

~~+ Now: Add sample dataset~~

+ Now: Add support for writing to anndata object
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
