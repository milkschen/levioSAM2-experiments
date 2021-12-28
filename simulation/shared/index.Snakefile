'''
Rules for indexing the source genome
'''
rule bt2_index_source:
    input:
        SOURCE_REF
    output:
        expand(os.path.join(
            DIR_IDX, f'{SOURCE_LABEL}' + '.{i}.bt2'),
            i = BT2_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, f'{SOURCE_LABEL}')
    threads: THREADS
    shell:
        '{BOWTIE2}-build --threads {threads} {input} {params}'

rule bwa_index_source:
    input:
        SOURCE_REF
    output:
        expand(os.path.join(
            DIR_IDX, f'{SOURCE_LABEL}' + '.{i}'),
            i = BWA_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, f'{SOURCE_LABEL}')
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
            DIR_IDX, f'{TARGET_LABEL}' + '.{i}.bt2'),
            i = BT2_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, f'{TARGET_LABEL}')
    threads: THREADS
    shell:
        '{BOWTIE2}-build --threads {threads} {input} {params}'

rule bwa_index_target:
    input:
        SOURCE_REF
    output:
        expand(os.path.join(
            DIR_IDX, f'{TARGET_LABEL}' + '.{i}'),
            i = BWA_IDX_ITEMS)
    params:
        os.path.join(DIR_IDX, f'{TARGET_LABEL}')
    shell:
        '{BWA} index -p {params} {input}'

rule check_index:
    input:
        expand(os.path.join(DIR_IDX, '{label}.{idx_item}.bt2'),
            idx_item = BT2_IDX_ITEMS,
            label = [SOURCE_LABEL, TARGET_LABEL]),
        expand(os.path.join(DIR_IDX, '{label}.{idx_item}'),
            idx_item = BWA_IDX_ITEMS,
            label = [SOURCE_LABEL, TARGET_LABEL])
    output:
        touch(temp(os.path.join(DIR, 'index.done')))
