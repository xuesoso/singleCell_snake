## Rules
"""
    Future goal: separate samples into groups of specified sizes
    Then we are going to perform STAR in each group while keeping the
    genome index in memory. This saves a lot of loading time.
    We don't want to sort by coordinate as it would increase htseq-count
    buffer requirement.
    We also keep all alignments with the best score, instead of keeping one
    of the best alignments.
"""

rule star:
    input:
        rules.star_index.output,
        get_all_fqgz
    output:
        "{all_samples}/Aligned.out.bam",
        "{all_samples}/Log.final.out"
    threads: 6
    params:
        name='star',
        partition=PART,
        mem=STAR_MEM,
        time=STAR_TIME
    shell:  "wdir=$(dirname {output[0]})/ && "
            "echo $wdir && "
            "STAR "
            "--genomeDir {input[0]} "
            "--readFilesIn {input[1]} {input[2]} "
            "--outSAMstrandField intronMotif "
            "--readFilesCommand gunzip -c "
            "--outFileNamePrefix $wdir "
            "--outSAMtype BAM Unsorted SortedByCoordinate "
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
