'''
Rules for one-pass alignment methods including grc, major and personalized.

grc and major are straightforward one-pass alignments, but the personalized
alignment method all are aligned to both hapA and hapB, and a merging step
is performed to generate the "best" outputs.

Checkpoint:
    temp(os.path.join(DIR, 'standard_onepass.done'))
'''
rule align_to_source:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_source.{i}.bt2'), i = IDX_ITEMS)
    params:
        DIR_IDX + 'chr{}_source'.format(CHROM)
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-source.bam'.format(CHROM))
    threads: THREADS
    shell:
        '{BOWTIE2} --threads {threads} -x {params} -1 {input.reads1} -2 {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output.sam}'

rule align_to_target:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_target.{i}.bt2'), i = IDX_ITEMS)
    params:
        DIR_IDX + 'chr{}_target'.format(CHROM)
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-target.bam'.format(CHROM))
    threads: THREADS
    shell:
        '{BOWTIE2} --threads {threads} -x {params} -1 {input.reads1} -2 {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output.sam}'

rule check_standard_onepass:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS, 'chr{}-source.bam'.format(CHROM)),
            INDIV = INDIV),
        expand(
            os.path.join(DIR_FIRST_PASS, 'chr{}-target.bam'.format(CHROM)),
            INDIV = INDIV),
        # samA = expand(
        #     os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA.sam'.format(CHROM)),
        #     INDIV = INDIV),
        # samB = expand(
        #     os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB.sam'.format(CHROM)),
        #     INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'alignment.done')))
