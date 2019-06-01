## Rules
rule make_transcript_annotation:
    input: REFERENCE_ANNOTATION
    output: TRANSCRIPT_ANNOTATION+'_transcript_annotation.tsv.gz'
    params:
        name='make_annotation',
        partition='quake,normal',
        mem='2000',
        time='20:00'
    run:
        make_annotation_dataframe(input[0], input[0])

rule make_gene_exp:
    input: "{outfile}/transcript_matrix/{outname}_merged_htseq.tab", rules.make_transcript_annotation.output
    output: "{outfile}/gene_matrix/{outname}_merged_htseq_gene.tab.gz"
    params:
        name='gene_expression',
        partition='quake,normal',
        mem='20000',
        time='30:00'
    run: convert_transcript_to_gene(input[0], input[1], output[0])
