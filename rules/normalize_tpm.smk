## Rules
rule normalize_tpm:
    input:
        htseq=rules.gzip_tables.output[0],
        star=rules.gzip_tables.output[1],
        ann=TRANSCRIPT_ANNOTATION+'_exon_annotation.tsv.gz',
    output:
        tpm="{outfile}/transcript_matrix/{outname}_merged_htseq.tpm.tab.gz"
    params:
        name='normalize_tpm',
        partition='quake,normal',
        mem='16000',
        time='20:00'
    run:
        tpm_normalize(input.htseq, input.ann, input.star, output.tpm)

rule make_tpm_exp:
    input:
        "{outfile}/transcript_matrix/{outname}_merged_htseq.tpm.tab.gz",
        rules.make_transcript_annotation.output
    output: "{outfile}/gene_matrix/{outname}_merged_htseq_gene.tpm.tab.gz"
    params:
        name='merge_tpm_expression',
        partition='quake,normal',
        mem='64000',
        time='30:00'
    run: convert_transcript_to_gene(input[0], input[1], output[0])
