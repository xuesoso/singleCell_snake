rule bam_to_vcf:
    """ Generates VCF for IGH region """
    input:
        bam='{base}/{sample}/STAR_output/Aligned.out.sorted.bam'
    output:
        '{base}/{sample}/STAR_output/full_variants.vcf.gz'
    params:
        name="vcf",
        partition=default_partition
    resources:
        mem_mb=5000
    conda:
        os.path.join(workflow.basedir, 'envs/miniconda.yaml')
    shell:
        "bcftools mpileup -Ou -f {REFERENCE_FASTA} {input.bam} |"
        " bcftools call -m -Ou |"
        " bcftools filter -e 'INFO/DP < 2' --output-type z -o {output}"
