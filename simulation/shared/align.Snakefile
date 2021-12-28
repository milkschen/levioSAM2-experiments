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
            DIR_IDX, f'{SOURCE_LABEL}' + '.{i}.bt2'), i = BT2_IDX_ITEMS)
    output:
        os.path.join(DIR_FIRST_PASS, f'bt2-{SOURCE_LABEL}.bam')
    params:
        os.path.join(DIR_IDX, f'{SOURCE_LABEL}')
    threads: THREADS
    shell:
        '{BOWTIE2} --threads {threads} -x {params} -1 {input.reads1} -2 {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule bwa_align_to_source:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, f'{SOURCE_LABEL}' + '.{i}'), i = BWA_IDX_ITEMS)
    output:
        os.path.join(DIR_FIRST_PASS, f'bwa-{SOURCE_LABEL}.bam')
    params:
        os.path.join(DIR_IDX, f'{SOURCE_LABEL}')
    threads: THREADS
    shell:
        '{BWA} mem -t {threads} {params} {input.reads1} {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule bt2_align_to_target:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, f'{TARGET_LABEL}' + '.{i}.bt2'), i = BT2_IDX_ITEMS)
    output:
        os.path.join(DIR_FIRST_PASS, f'bt2-{TARGET_LABEL}.bam')
    params:
        os.path.join(DIR_IDX, f'{TARGET_LABEL}')
    threads: THREADS
    shell:
        '{BOWTIE2} --threads {threads} -x {params} -1 {input.reads1} -2 {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule bwa_align_to_target:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq',
        idx = expand(os.path.join(
            DIR_IDX, f'{TARGET_LABEL}' + '.{i}'), i = BWA_IDX_ITEMS)
    output:
        os.path.join(DIR_FIRST_PASS, f'bwa-{TARGET_LABEL}.bam')
    params:
        os.path.join(DIR_IDX, f'{TARGET_LABEL}')
    threads: THREADS
    shell:
        '{BWA} mem -t {threads} {params} {input.reads1} {input.reads2} | '
        '{SAMTOOLS} view -hb -o {output}'

rule check_alignment:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS, '{aln}-{label}.bam'),
            INDIV = INDIV, label = [SOURCE_LABEL, TARGET_LABEL],
            aln = ['bwa', 'bt2']),
    output:
        touch(temp(os.path.join(DIR, 'alignment.done')))
