""" For TPM normalization, we need to calculate the effective insert/fragment size """
rule run_picard:
    input:
        bam=rules.star.output[0]
    output:
        insert_out='{all_samples}/insert.out'
    params:
        name='average_insert_length',
        picard=PICARD,
        partition=PART,
        mem='8000',
        time='00:15:00'
    shell:
        """
        touch {output.insert_out}
        {params.picard} CollectInsertSizeMetrics -I {input.bam} -H {output.insert_out}.hist -O {output.insert_out}
        """

rule merge_insert_size:
    input:
        insert_out=rules.run_picard.output.insert_out,
        STAR=rules.star.output[1]
    output:
        STAR='{all_samples}/STAR.out'
    params:
        name='average_insert_length',
        partition=PART,
        mem='4000',
        time='00:4:00'
    shell:
        """
        cp {input.STAR} {output.STAR}
        if [ -s {input.insert_out} ]
        then
            echo -n 'average insert length | ' >> {output.STAR} && grep -A 1 'MEAN_INSERT_SIZE' {input.insert_out} | sed -n '2p' | cut -f6 | tr -d '\n' >> {output.STAR}
            echo -n '\nmedian insert length | ' >> {output.STAR} && grep -A 1 'MEDIAN_INSERT_SIZE' {input.insert_out} | sed -n '2p' | cut -f1 | tr -d '\n' >> {output.STAR}
        else
            echo 'average insert length | 0' >> {output.STAR}
            echo 'median insert length | 0' >> {output.STAR}
        fi
        """
