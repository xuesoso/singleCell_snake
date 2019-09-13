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
    STAR_TIME = '8:00:00'
else:
    SPARSITY = 2
    STAR_MEM = 30000
    STAR_TIME = '8:00:00'

## Cluster parameters
PART = config['partition']
# CHUNKSIZE = config['chunksize']

## STAR keep one-best aligned or keep all best aligned
if "star_keep" not in config.keys():
    STAR_KEEP = 'one'
else:
    STAR_KEEP = config['star_keep']

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

## Define rules
include: "./rules/genome_index.smk"

""" Align reads using STAR. """
if STAR_KEEP == 'one':
    """ Keep only one of the best aligned reads as primary.
    Mark other reads as secondary. In Htseq, these secondary reads won't be
    counted. This is deterministic.
    """
    include: './rules/STAR_keep_best_one.smk'

elif STAR_KEEP == 'all':
    """ Mark all of best aligned reads as primary reads. They can be all be
    counted by Htseq as contributing to one over N to all N multimapped
    regions.
    """
    include: './rules/STAR_keep_all_best.smk'

elif STAR_KEEP == 'unmapped':
    """ STAR align the Unmapped read pairs.
    """
    include: './rules/STAR_unmapped.smk'

else:
    raise ValueError("STAR_KEEP must be either 'one' or 'all'")


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

if STAR_KEEP != 'unmapped':
    rule all:
        input:
            expand("{outfile}/gene_matrix/{outname}_merged_htseq_gene.tab.gz", outfile=outfile, outname=outname),
            expand("{outfile}/star_matrix/{outname}_merged_star.tab.gz", outfile=outfile, outname=outname),
            expand("{outfile}/snp_matrix/{outname}_merged_vcf.tab.gz", outfile=outfile, outname=outname)
        params:
            name='all',
            partition='quake,normal',
            mem='1024',
            time='1:00'
else:
    rule all:
        input: expand("{all_samples}/Unmapped.Log.final.out", all_samples=all_samples)
        params:
            name='all',
            partition='quake,normal',
            mem='1024',
            time='1:00'
