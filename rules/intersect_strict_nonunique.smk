""" Use intersection-strict mode """
""" We count all multimapping reads by adding a fraction of the counts """
rule htseq:
    input:
        bam=rules.star.output[0]
    output: '{all_samples}/htseq.tab'
    params:
        name='htseq',
        partition=PART,
        mem='5000',
        time='01:00:00'
    shell: "htseq-count -s no -f bam -m intersection-strict "
            "-i Parent "
            "--nonunique all "
            "--stranded no "
            "-r name "
            "--secondary-alignments ignore "
            "--supplementary-alignments ignore "
            "{input.bam} {REFERENCE_ANNOTATION} > {output}"
