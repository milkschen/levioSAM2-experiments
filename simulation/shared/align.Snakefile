'''
Rules for alignment jobs.

Checkpoint:
    temp(os.path.join(DIR, 'alignment.done'))
'''
rule bt2_align_to_source:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_source.{i}.bt2'), i = BT2_IDX_ITEMS)
    params:
        DIR_IDX + 'chr{}_source'.format(CHROM)
    output:
        os.path.join(DIR_FIRST_PASS, 'bt2-chr{}-source.bam'.format(CHROM))
    threads: THREADS
    shell:
        '{BOWTIE2} --threads {threads} -x {params} -1 {input.reads1} -2 {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule bwa_align_to_source:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_source.{i}'), i = BWA_IDX_ITEMS)
    params:
        DIR_IDX + 'chr{}_source'.format(CHROM)
    output:
        os.path.join(DIR_FIRST_PASS, 'bwa-chr{}-source.bam'.format(CHROM))
    threads: THREADS
    shell:
        '{BWA} mem -t {threads} {params} {input.reads1} {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule bt2_align_to_target:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_target.{i}.bt2'), i = BT2_IDX_ITEMS)
    params:
        DIR_IDX + 'chr{}_target'.format(CHROM)
    output:
        os.path.join(DIR_FIRST_PASS, 'bt2-chr{}-target.bam'.format(CHROM))
    threads: THREADS
    shell:
        '{BOWTIE2} --threads {threads} -x {params} -1 {input.reads1} -2 {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule bwa_align_to_target:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_target.{i}'), i = BWA_IDX_ITEMS)
    params:
        DIR_IDX + 'chr{}_target'.format(CHROM)
    output:
        os.path.join(DIR_FIRST_PASS, 'bwa-chr{}-target.bam'.format(CHROM))
    threads: THREADS
    shell:
        '{BWA} mem -t {threads} {params} {input.reads1} {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule check_standard_onepass:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS, 'bt2-chr{}-source.bam'.format(CHROM)),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, 'bt2-chr{}-target.bam'.format(CHROM)),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, 'bwa-chr{}-source.bam'.format(CHROM)),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, 'bwa-chr{}-target.bam'.format(CHROM)),
            INDIV = INDIV),
        # samA = expand(
        #     os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA.sam'.format(CHROM)),
        #     INDIV = INDIV),
        # samB = expand(
        #     os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB.sam'.format(CHROM)),
        #     INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'alignment.done')))
