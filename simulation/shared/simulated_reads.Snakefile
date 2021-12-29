rule extract_genotypes:
    input:
        vcf = PHASED_VCF
    output:
        vcf = os.path.join(DIR, 'simulation/{s}/' + 'chr{}-per'.format(CHROM) + '.vcf.gz')
    shell:
        '{BCFTOOLS} view -s {wildcards.s} -v snps {input.vcf} | '
        '{BCFTOOLS} norm -d all -O z -o {output.vcf}'

rule simulate_reads:
    input:
        vcf = PREFIX_PER + '.vcf.gz',
        ref = SIMULATION_REF
    output:
        fq1 = temp(PREFIX_PER + '_1.fq'),
        fq2 = temp(PREFIX_PER + '_2.fq'),
        sam = temp(PREFIX_PER + '.sam')
    params:
        num = NUM_SIM_READS
    threads: MAX_SYSTEM_THREADS
    shell:
        '{MASON2} --num-threads {threads} -ir {input.ref} -n {params.num} '
        '-o {output.fq1} -or {output.fq2} -oa {output.sam} -iv {input.vcf}'

rule simulation_compress_fq:
    input: PREFIX_PER + '_{suffix}.fq'
    output: PREFIX_PER + '_{suffix}.fq.gz'
    shell: 'gzip {input}'

rule simulation_compress_sam:
    input: PREFIX_PER + '.sam'
    output: PREFIX_PER + '.bam'
    shell: '{SAMTOOLS} view -hb -o {output} {input}'


rule check_simulation:
    input:
        expand(
            PREFIX_PER + '_{seg}.fq.gz',
            INDIV = INDIV, seg = ['1', '2']),
        expand(
            PREFIX_PER + '.bam',
            INDIV = INDIV)
    output:
        touch(temp(os.path.join(DIR, 'simulation.done')))

