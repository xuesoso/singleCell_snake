## Rules
""" Create a genome reference index """
rule star_index:
    input:
        fasta=REFERENCE_FASTA,
        annotation=REFERENCE_ANNOTATION
    output:
        ### Since Snakemake 5.2 we need to specify directory output
        directory(REFERENCE_INDEX)
    params:
        name='star_index',
        partition=PART,
        mem=64000,
        time='1:00:00',
        threads=16,
        sparsity=SPARSITY
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate --runThreadN {params.threads} "
        "--genomeSAsparseD {params.sparsity} --genomeFastaFiles {input.fasta} "
        "--genomeDir {output} --sjdbGTFfile {input.annotation} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"
