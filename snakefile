include: "Snakefile_utils.py"
import os

## Input locations
REFERENCE_ANNOTATION = config['annotation']
TRANSCRIPT_ANNOTATION = transcript_annotation_name(REFERENCE_ANNOTATION)
# REFERENCE_INDEX = config['refIndex']
REFERENCE_FASTA = config['refFasta']
REFERENCE_INDEX = '.'.join(REFERENCE_FASTA.split('.')[:-1]) + '.starIndex'
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
print(all_samples)

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
include: "./rules/snp.smk"
include: "./rules/genome_index.smk"

rule all:
    input:
        expand("{outfile}/gene_matrix/{outname}_merged_htseq_gene.tab.gz", outfile=outfile, outname=outname),
        expand("{outfile}/star_matrix/{outname}_merged_star.tab.gz", outfile=outfile, outname=outname),
        # rules.bam_to_vcf.output
        expand("{all_samples}/full_variants.vcf.gz", all_samples=all_samples)
    params:
        name='all',
        partition='quake,normal',
        mem='1024',
        time='1:00'
