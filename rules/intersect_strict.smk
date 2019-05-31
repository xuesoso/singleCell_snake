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

