## Rules
""" Create a genome reference index """
rule star_index:
    input:
        REFERENCE_FASTA,
        REFERENCE_ANNOTATION
    output:
        REFERENCE_INDEX
    params:
        name='star_index',
        partition=PART,
        mem=32000,
        time='1:00:00',
        threads=16,
        sparsity=2
    shell:
        "mkdir -p {output} && "
        "STAR --runMode genomeGenerate --runThreadN {params.threads} "
        "--genomeSAsparseD {params.sparsity} --genomeFastaFiles {input[0]} "
        "--genomeDir {output} --sjdbGTFfile {input[1]} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"
