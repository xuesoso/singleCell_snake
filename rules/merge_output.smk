## Rules
rule merge_htseq:
    input: expand("{all_samples}/htseq.tab", all_samples=all_samples)
    output: "{outfile}/transcript_matrix/{outname}_merged_htseq.tab"
    params:
        name='merge_htseq',
        partition='quake,normal',
        mem='30000',
        time='1:00:00'
    run: merge_htseq_tables(input, output[0])

rule merge_star:
    input: expand("{all_samples}/Log.final.out", all_samples=all_samples)
    output: "{outfile}/star_matrix/{outname}_merged_star.tab"
    params:
        name='merge_star',
        partition='quake,normal',
        mem='30000',
        time='1:00:00'
    run: merge_star_tables(input, output[0])

rule gzip_tables:
    input:
        rules.merge_htseq.output,
        rules.merge_star.output
    output:
        "{outfile}/transcript_matrix/{outname}_merged_htseq.tab.gz",
        "{outfile}/star_matrix/{outname}_merged_star.tab.gz"
    params:
        name='gzip_tables',
        partition='quake,normal',
        mem='10000',
        time='10:00'
    shell:
        "gzip {input}"
