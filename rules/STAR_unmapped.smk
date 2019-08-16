## Rules
"""
    We will gather the mated pairs of Unmapped.out.mate1 and Unmapped.out.mate2, 
    usually produced by a previous STAR run. We can then align these Unmapped reads
    to a genome reference index containing TSO / PCR primers to check self-homology.
"""

rule star:
    input:
        rules.star_index.output,
        get_all_fqgz
    output:
        "{all_samples}/Unmapped.Aligned.out.bam",
        "{all_samples}/Unmapped.Log.final.out",
        "{all_samples}/Unmapped.Aligned.sortedByCoord.out.bam"
    threads: 6
    params:
        name='star_unmapped',
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
            "--outFilterType BySJout "
            "--outFilterMultimapNmax 20 "
            "--outSAMprimaryFlag OneBestScore "
            "--outFilterScoreMinOverLread 0.4 "
            "--outFilterMatchNminOverLread 0.4 "
            "--outFilterMismatchNmax 999 "
            "--outFilterMismatchNoverLmax 0.04 "
            "--alignIntronMin 20 "
            "--alignIntronMax 1000000 "
            "--alignMatesGapMax 1000000 "
            "--alignSJoverhangMin 8 "
            "--alignSJDBoverhangMin 1"
