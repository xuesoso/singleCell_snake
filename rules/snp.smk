"""
Scripts originally provided by Derek Croote
Tested to work on current pipeline 6/26/2019
"""

rule bam_to_vcf:
    """ Generates variable snp calling """
    input:
        fasta=REFERENCE_FASTA,
        bam=rules.star.output[0]
    output:
        "{all_samples}/full_variants.vcf.gz"
    params:
        name="vcf",
        partition=PART,
        mem=9000,
        time='1:00:00'
    shell:
        "bcftools mpileup -Ou -f {input.fasta} {input.bam} |"
        " bcftools call -mv -Ou |"
        # " bcftools call -m -Ou |" ## output all sites
        " bcftools filter -e 'INFO/DP < 2' --output-type z -o {output}"
