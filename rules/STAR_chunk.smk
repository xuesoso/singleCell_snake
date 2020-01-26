## Rules
""" Future goal: separate samples into groups of specified sizes """
""" Then we are going to perform STAR in each group while keeping the """
""" genome index in memory. This saves a lot of loading time """
rule create_chunks:
    input:
        get_all_fqgz
    output:
        "{chunks}.in"
    run:
        import os
        samples_to_reads = {}
        for read_pairs in input:
            sample = os.path.dirname(read_pairs[0])
            samples_to_reads[sample] = read_pairs
        new_dict = {}
        for k in chunks_to_samples:
            new_dict[k] = samples_to_reads[chunks_to_samples[k]]
        write_chunks(new_dict, postfix='.in')

rule star:
    input:
        rules.star_index.output,
        "{chunks}.in"
        # expand("{chunks}.in", chunks=chunks)
        # get_all_fqgz
    output:
        # "{all_samples}/Aligned.sortedByCoord.out.bam",
        # "{all_samples}/Log.final.out"
        # temp("{chunks}")
        "{chunks}.done"
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
            "--twopassMode Basic"
            # "--alignIntronMin 20 "
            # "--alignIntronMax 1000000 "
            # "--alignMatesGapMax 1000000 "
            # "--alignSJoverhangMin 8 "
            # "--alignSJDBoverhangMin 1"
