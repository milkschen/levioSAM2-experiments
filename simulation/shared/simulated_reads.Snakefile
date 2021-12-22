'''
Rules for simulate reads from personalized genomes
'''
rule simulate_reads:
    input:
        hapA = PREFIX_PER + '_hapA.fa',
        hapB = PREFIX_PER + '_hapB.fa'
    output:
        readsA1 = temp(PREFIX_PER + '_hapA_1.fq'),
        readsA2 = temp(PREFIX_PER + '_hapA_2.fq'),
        readsB1 = temp(PREFIX_PER + '_hapB_1.fq'),
        readsB2 = temp(PREFIX_PER + '_hapB_2.fq'),
        samA = PREFIX_PER + '_hapA.sam',
        samB = PREFIX_PER + '_hapB.sam'
    params:
        num = NUM_SIM_READS,
        prefix = PREFIX_PER,
    #: set to MAX_SYSTEM_THREADS to avoid errors due to shared temp files
    threads: MAX_SYSTEM_THREADS
    shell:
        '{MASON2} --num-threads {threads} -ir {input.hapA} -n {params.num} '
        '-o {params.prefix}_hapA_1.fq -or {params.prefix}_hapA_2.fq '
        '-oa {params.prefix}_hapA.sam --read-name-prefix "{params.prefix}_hapA_simulated.";'
        '{MASON2} --num-threads {threads} -ir {input.hapB} -n {params.num} '
        '-o {params.prefix}_hapB_1.fq -or {params.prefix}_hapB_2.fq '
        '-oa {params.prefix}_hapB.sam --read-name-prefix "{params.prefix}_hapB_simulated.";'

rule merge_simulated_reads:
    input:
        readsA1 = PREFIX_PER + '_hapA_1.fq',
        readsA2 = PREFIX_PER + '_hapA_2.fq',
        readsB1 = PREFIX_PER + '_hapB_1.fq',
        readsB2 = PREFIX_PER + '_hapB_2.fq',
    output:
        reads1 = PREFIX_PER + '_1.fq',
        reads2 = PREFIX_PER + '_2.fq'
    shell:
        'cat {input.readsA1} > {output.reads1};'
        'cat {input.readsA2} > {output.reads2};'
        'cat {input.readsB1} >> {output.reads1};'
        'cat {input.readsB2} >> {output.reads2};'

rule merge_simulated_sam:
    input:
        samA = PREFIX_PER + '_hapA.sam',
        samB = PREFIX_PER + '_hapB.sam'
    output:
        sam = temp(PREFIX_PER + '.sam'),
        sam1 = PREFIX_PER + '_1.sam',
        sam2 = PREFIX_PER + '_2.sam'
    shell:
        'grep ^@ {input.samA} > {output.sam};'
        #: get headers from B
        'grep ^@ {input.samB} | tail -n +2 >> {output.sam};'
        'grep -v ^@ {input.samA} >> {output.sam};'
        'grep -v ^@ {input.samB} >> {output.sam};'
        '{SAMTOOLS} view -h -f 64 {output.sam} -o {output.sam1};'
        '{SAMTOOLS} view -h -f 128 {output.sam} -o {output.sam2}'

rule check_simulation:
    input:
        expand(
            PREFIX_PER + '_{seg}.{type}',
            INDIV = INDIV,
            seg = ['1', '2'],
            type = ['fq', 'sam']
        )
    output:
        touch(temp(os.path.join(DIR, 'simulation.done')))

