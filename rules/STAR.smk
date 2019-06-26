## Rules
""" Future goal: separate samples into groups of specified sizes """
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
        time='12:00:00'
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
