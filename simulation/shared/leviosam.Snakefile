rule bt2_lift_source_to_target:
    input:
        bam = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-source.bam'),
        clft = CLFT,
        tref = TARGET_REF
    output:
        commit = temp(os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-committed.bam')),
        defer = temp(os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-deferred.bam')),
        unliftable = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-unliftable.bam')
    params:
        prefix = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}'),
        options = LEVIOSAM_OPTIONS,
        g = MAX_CHAIN_GAP
    threads: THREADS
    shell:
        '{LEVIOSAM} lift -C {input.clft} -a {input.bam} -t {threads} -p {params.prefix} '
        '-m -f {input.tref} -O bam -G {params.g} {params.options}'

rule bt2_collate:
    input:
        commit = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-committed.bam'),
        defer = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-deferred.bam'),
    output:
        commit = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-committed.bam'),
        defer = temp(os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred.bam')),
        fq1 = temp(os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-R1.fq.gz')),
        fq2 = temp(os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-R2.fq.gz')),
    params:
        prefix = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired'),
    shell:
        '{LEVIOSAM} collate -a {input.commit} -b {input.defer} -p {params.prefix}'

rule bt2_sort_defer:
    input:
        bam = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred.bam'),
    output:
        bam = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-sorted_n.bam'),
    threads: int(THREADS * 0.2)
    shell:
        '{SAMTOOLS} sort -@ {threads} -n -o {output.bam} {input.bam}'

rule bt2_realign:
    input:
        fq1 = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-R1.fq.gz'),
        fq2 = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-R2.fq.gz'),
        idx = expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_target.{i}.bt2'), i = BT2_IDX_ITEMS)
    output:
        bam = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-realigned-sorted_n.bam')
    params:
        idx = DIR_IDX + 'chr{}_target'.format(CHROM),
        rg = BT2_READ_GROUP,
        t_80 = int(THREADS * 0.8),
        t_20 = THREADS - int(THREADS * 0.8)
    threads: THREADS
    shell:
        '{BOWTIE2} {params.rg} -p {params.t_80} -x {params.idx} '
        '-1 {input.fq1} -2 {input.fq2} | '
        '{SAMTOOLS} sort -@ {params.t_20} -n -o {output.bam}'

rule bt2_cherry_pick:
    input:
        defer = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-sorted_n.bam'),
        realn = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-realigned-sorted_n.bam')
    output:
        bam = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-cherry_picked.bam')
    shell:
        '{LEVIOSAM} cherry_pick -s source:{input.defer} -s target:{input.realn} -m -o {output.bam}'

rule bt2_merge_and_sort:
    input:
        commit = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-committed.bam'),
        cherry_picked = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-paired-deferred-cherry_picked.bam')
    output:
        bam = os.path.join(DIR_FIRST_PASS, f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-final.bam')
    threads: int(THREADS * 0.4)
    shell:
        '{SAMTOOLS} cat {input.commit} {input.cherry_picked} | '
        '{SAMTOOLS} sort -@ {threads} -o {output.bam}'

rule check_leviosam:
    input:
        expand(
            os.path.join(DIR_FIRST_PASS,
            f'bt2-chr{CHROM}-{SOURCE_LABEL}_to_{TARGET_LABEL}-final.bam'),
            INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'leviosam.done')))
