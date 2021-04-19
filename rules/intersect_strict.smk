""" Use intersection-strict mode """
rule htseq:
    input:
        bam=rules.star.output[0]
    output: '{all_samples}/htseq.tab'
    params:
        name='htseq',
        partition=PART,
        mem='8000',
        time='01:00:00'
    shell: "python -m HTSeq.scripts.count -s no -f bam -m intersection-strict "
            "-r name "
            "-i Parent "
            "{input.bam} "
            "{REFERENCE_ANNOTATION} > {output}"
