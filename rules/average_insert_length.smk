""" For TPM normalization, we need to calculate the effective insert/fragment size """
rule average_insert_length:
    input:
        bam=rules.star.output[0],
        STAR=rules.star.output[1],
    output:
        STAR='{all_samples}/STAR.out',
        insert_hist='{all_samples}/insert.hist.out',
        insert_out='{all_samples}/insert.out'
    params:
        name='average_insert_length',
        picard=PICARD,
        partition=PART,
        mem='8000',
        time='00:15:00'
    shell:
        """
        cp {input.STAR} {output.STAR}
        echo -n 'average insert length | ' >> {output.STAR}
        {params.picard} CollectInsertSizeMetrics -I {input.bam} -H {output.insert_hist} -O {output.insert_out} && grep -A 1 'MEAN_INSERT_SIZE' {output.insert_out} | sed -n '2p' | cut -f6 | tr -d '\n' >> {output.STAR}
        """
