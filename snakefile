include: "Snakefile_utils.py"
import os

## Input locations
REFERENCE_ANNOTATION = config['annotation']
TRANSCRIPT_ANNOTATION = transcript_annotation_name(REFERENCE_ANNOTATION)
REFERENCE_INDEX = config['refIndex']
INPUTFILE = config['inputDir']
PLATES = config['plates']

## Rule parameters
HTSEQ_MODE = config['htseq_mode']

## Cluster parameters
PART = config['partition']
CHUNKSIZE = config['chunksize']

## Define wildcards for inputs
root_dir = INPUTFILE
plate = PLATES
all_samples = []
for f in PLATES:
    base = os.path.join(root_dir, f)
    for j in (os.listdir(base)):
        all_samples.extend([os.path.join(base, j)])

## Define wildcards for outputs
outfile = config['outputDir']
if 'outname' not in config.keys():
    outname = "final"
else:
    outname = config['outname']

include: './rules/STAR.smk'

""" Count reads mapping to features using htseq """
if HTSEQ_MODE == 'union':
    """ Use union mode """
    include: './rules/union.smk'

elif HTSEQ_MODE == 'intersect_strict':
    """ Use intersection-strict mode """
    include: './rules/intersect_strict.smk'

else:
    raise ValueError("HTSEQ_MODE must be either 'union' or 'intersect_strict'")

include: "./rules/merge_output.smk"
include: "./rules/feature_to_gene.smk"

rule all:
    input:
        expand("{outfile}/gene_matrix/{outname}_merged_htseq_gene.tab.gz", outfile=outfile, outname=outname),
        rules.merge_star.output
    params: name='all', partition='quake,normal', mem='1024'
