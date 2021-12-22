'''
Rules for indexing the source genome
'''
rule bt2_index_source:
    input:
        SOURCE_REF
    output:
        expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_source.{i}.bt2'),
            i = BT2_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, 'chr{}_source'.format(CHROM))
    threads: THREADS
    shell:
        '{BOWTIE2}-build --threads {threads} {input} {params}'

rule bwa_index_source:
    input:
        SOURCE_REF
    output:
        expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_source.{i}'),
            i = BWA_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, 'chr{}_source'.format(CHROM))
    threads: THREADS
    shell:
        '{BWA} index -p {params} {input}'

'''
Rules for indexing the target genome
'''
rule bt2_index_target:
    input:
        TARGET_REF
    output:
        expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_target.{i}.bt2'),
            i = BT2_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, 'chr{}_target'.format(CHROM))
    threads: THREADS
    shell:
        '{BOWTIE2}-build --threads {threads} {input} {params}'

rule bwa_index_target:
    input:
        SOURCE_REF
    output:
        expand(os.path.join(
            DIR_IDX, 'chr{}'.format(CHROM) + '_target.{i}'),
            i = BWA_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, 'chr{}_target'.format(CHROM))
    threads: THREADS
    shell:
        '{BWA} index -p {params} {input}'

rule check_index:
    input:
        expand(os.path.join(DIR_IDX, 'chr{}'.format(CHROM) + '_{label}.{idx_item}.bt2'),
            idx_item = BT2_IDX_ITEMS,
            label = ['source', 'target']),
        expand(os.path.join(DIR_IDX, 'chr{}'.format(CHROM) + '_{label}.{idx_item}'),
            idx_item = BWA_IDX_ITEMS,
            label = ['source', 'target'])
    output:
        touch(temp(os.path.join(DIR, 'index.done')))
