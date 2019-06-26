"""
Scripts originally provided by Derek Croote
Not tested on current rule sets.
"""

rule bam_to_vcf:
    """ Generates variable snp calling """
    input:
        # bam='{base}/{sample}/STAR_output/Aligned.out.sorted.bam'
        bam=rules.star.output[0]
    output:
        "{all_samples}/full_variants.vcf.gz"
    params:
        name="vcf",
        partition=PART,
        mem=9000,
        time='1:00:00'
    shell:
        "bcftools mpileup -Ou -f {REFERENCE_FASTA} {input.bam} |"
        " bcftools call -m -Ou |"
        " bcftools filter -e 'INFO/DP < 2' --output-type z -o {output}"
