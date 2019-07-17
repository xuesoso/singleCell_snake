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
    shell: "samtools view {input.bam} | htseq-count -s no -f sam -m "
            "intersection-strict -r name -i Parent - {REFERENCE_ANNOTATION} > "
            "{output}"
