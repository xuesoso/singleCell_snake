""" Use intersection-strict mode """
""" We count all multimapping reads by adding a fraction of the counts """
rule htseq:
    input:
        bam=rules.star.output[0]
    output: '{all_samples}/htseq.tab'
    params:
        name='htseq',
        partition=PART,
        mem='8000',
        time='01:00:00'
    shell: "samtools view {input.bam} | "
            "htseq-count -s no -f sam -m intersection-strict "
            "-i Parent "
            "--nonunique all "
            "--stranded no "
            "-r name "
            "--secondary-alignments ignore "
            "--supplementary-alignments ignore "
            "- {REFERENCE_ANNOTATION} > {output}"
