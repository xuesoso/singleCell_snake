include: "Snakefile_utils.py"
import os

## Input locations
REFERENCE_ANNOTATION = config['annotation']
TRANSCRIPT_ANNOTATION = transcript_annotation_name(REFERENCE_ANNOTATION)
REFERENCE_FASTA = config['refFasta']
REFERENCE_INDEX = '.'.join(REFERENCE_FASTA.split('.')[:-1]) + '.starIndex'
INPUTFILE = config['inputDir']
PLATES = config['plates']

## Rule parameters
HTSEQ_MODE = config['htseq_mode']

## Set sparsity of genome index (2 if more than 1G, 1 if less than 1G)
size_of_fasta = os.stat(REFERENCE_FASTA).st_size <= 107374182
if size_of_fasta:
    SPARSITY = 1
    STAR_MEM = 8000
    STAR_TIME = '12:00:00'
else:
    SPARSITY = 2
    STAR_MEM = 30000
    STAR_TIME = '12:00:00'

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

include: "./rules/genome_index.smk"
include: './rules/STAR.smk'

""" Count reads mapping to features using htseq """
if HTSEQ_MODE == 'unique':
    """ Use intersection-strict mode. Count only uniquely mapped. """
    include: './rules/intersect_strict.smk'

elif HTSEQ_MODE == 'multimap':
    """ Use intersection-strict mode. Count all equally well-mapped.
    Add a fraction of the features aligned """
    include: './rules/intersect_strict_nonunique.smk'

else:
    raise ValueError("HTSEQ_MODE must be either 'unique' or 'multimap'")

include: "./rules/merge_output.smk"
include: "./rules/feature_to_gene.smk"
include: "./rules/snp.smk"

rule all:
    input:
        expand("{outfile}/gene_matrix/{outname}_merged_htseq_gene.tab.gz", outfile=outfile, outname=outname),
        expand("{outfile}/star_matrix/{outname}_merged_star.tab.gz", outfile=outfile, outname=outname),
        expand("{all_samples}/variants.vcf.gz", all_samples=all_samples)
    params:
        name='all',
        partition='quake,normal',
        mem='1024',
        time='1:00'
