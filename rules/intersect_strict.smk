""" Use intersection-strict mode """
rule htseq:
    # input:  '{dir}/Aligned.out.bam'
    input:
        bam=rules.star_chunks.output
    output: '{dir}/htseq.tab'
    params:
        name='htseq',
        partition=PART,
        mem='16000',
        time='30:00:00'
    shell: "htseq-count -s no -r pos -f bam -m intersection-strict "
            "-i Parent {input.bam} {REFERENCE_ANNOTATION} > {output}"
