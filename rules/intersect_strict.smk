""" Use intersection-strict mode """
rule htseq:
    input:
        bam=rules.star.output
    output: '{all_samples}/htseq.tab'
    params:
        name='htseq',
        partition=PART,
        mem='10000',
        time='01:00:00'
    shell: "htseq-count -s no -f bam -m intersection-strict -r name "
            "-i Parent {input.bam} {REFERENCE_ANNOTATION} > {output}"

