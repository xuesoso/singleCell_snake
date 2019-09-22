"""
Scripts originally provided by Derek Croote
Tested to work on current pipeline 6/26/2019
"""

rule bam_to_vcf:
    """ Generates variable snp calling """
    input:
        fasta=REFERENCE_FASTA,
        sorted_bam=rules.star.output[2]
    output:
        "{all_samples}/variants.vcf.gz",
        "{all_samples}/variants.vcf.gz.csi"
    params:
        name="vcf",
        partition=PART,
        mem=9000,
        time='1:00:00'
    shell:
        "bcftools mpileup -Ou -f {input.fasta} {input.sorted_bam} |"
        " bcftools call -mv -Ou |"
        # " bcftools call -m -Ou |" ## output all sites
        " bcftools filter -e 'INFO/DP < 2' --output-type z -o {output[0]}"
        " && bcftools index -f {output[0]}"
