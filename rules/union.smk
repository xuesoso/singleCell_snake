""" Use union mode """
rule htseq:
    input:
        bam=rules.star.output
    output: '{all_samples}/htseq.tab'
    params:
        name='htseq',
        partition=PART,
        mem='5000',
        time='01:00:00'
    shell: "htseq-count -s no -f bam -m union "
            "-i Parent "
            "--nonunique all "
            "--stranded no "
            "-r name "
            "--secondary-alignments ignore "
            "--supplementary-alignments ignore "
            "{input.bam} {REFERENCE_ANNOTATION} > {output}"
