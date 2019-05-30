include: "Snakefile_utils.py"
import os

## Input locations
REFERENCE_ANNOTATION = config['annotation']
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
# groups, chunkid = chunks(all_samples, CHUNKSIZE)

## Define wildcards for outputs
outfile = config['outputDir']
if 'outname' not in config.keys():
    outname = "final"
else:
    outname = config['outname']

## Rules
""" We are going to separate samples into groups of specified sizes """
""" Then we are going to perform STAR in each group while keeping the """
""" genome index in memory. This saves a lot of loading time """
rule star:
    input:
        REFERENCE_INDEX,
        get_all_fqgz
    output:
        "{all_samples}/Aligned.sortedByCoord.out.bam",
        "{all_samples}/Log.final.out"
    threads: 6
    params:
        name='star',
        partition=PART,
        mem=30000,
    shell:  "wdir=$(dirname {output[0]})/ && "
            "echo $wdir && "
            "STAR "
            "--genomeDir {input[0]} "
            "--readFilesIn {input[1]} {input[2]} "
            "--outSAMstrandField intronMotif "
            "--readFilesCommand gunzip -c "
            "--outFileNamePrefix $wdir "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes NH HI AS NM MD "
            "--outReadsUnmapped Fastx "
            "--clip3pAdapterSeq CTGTCTCTTATACACATCT "
            "--outFilterType BySJout "
            "--outFilterMultimapNmax 20 "
            "--outSAMprimaryFlag AllBestScore "
            "--outFilterScoreMinOverLread 0.4 "
            "--outFilterMatchNminOverLread 0.4 "
            "--outFilterMismatchNmax 999 "
            "--outFilterMismatchNoverLmax 0.04 "
            "--alignIntronMin 20 "
            "--alignIntronMax 1000000 "
            "--alignMatesGapMax 1000000 "
            "--alignSJoverhangMin 8 "
            "--twopassMode Basic "
            "--alignSJDBoverhangMin 1"

""" Count reads mapping to features using htseq """
if HTSEQ_MODE == 'union':
    """ Use union mode """
    rule htseq:
        input:
            bam=rules.star.output
        output: '{all_samples}/htseq.tab'
        params:
            name='htseq',
            partition=PART,
            mem='5000',
        shell: "htseq-count -s no -r pos -f bam -m intersection-strict "
                "-i Parent {input.bam} {REFERENCE_ANNOTATION} > {output}"

elif HTSEQ_MODE == 'intersect_strict':
    """ Use intersection-strict mode """
    rule htseq:
        input:
            bam=rules.star.output
        output: '{all_samples}/htseq.tab'
        params:
            name='htseq',
            partition=PART,
            mem='5000',
        shell: "htseq-count -s no -r pos -f bam -m intersection-strict "
                "-i Parent {input.bam} {REFERENCE_ANNOTATION} > {output}"
else:
    raise ValueError("HTSEQ_MODE must be either 'union' or 'intersect_strict'")

rule merge_htseq:
    input: expand("{all_samples}/htseq.tab", all_samples=all_samples)
    output: expand("{outfile}/transcript_matrix/{outname}_merged_htseq.tab", outfile=outfile, outname=outname)
    params: name='merge_htseq', partition='quake,normal', mem='30000', time='1:00:00'
    run: merge_htseq_tables(input, str(output[0]))

rule merge_star:
    input: expand("{all_samples}/Log.final.out", all_samples=all_samples)
    output: expand("{outfile}/star_matrix/{outname}_merged_star.tab", outfile=outfile, outname=outname)
    params: name='merge_star', partition='quake,normal', mem='30000', time='1:00:00'
    run: merge_star_tables(input, str(output[0]))

rule all:
    input:
        rules.merge_htseq.output, rules.merge_star.output
    params: name='all', partition='quake,normal', mem='1024'
