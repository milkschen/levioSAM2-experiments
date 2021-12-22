'''
Rules for one-pass alignment methods including grc, major and personalized.

grc and major are straightforward one-pass alignments, but the personalized
alignment method all are aligned to both hapA and hapB, and a merging step
is performed to generate the "best" outputs.

Checkpoint:
    temp(os.path.join(DIR, 'standard_onepass.done'))
'''
rule align_to_grc:
    input:
        reads1 = PREFIX_PER + '_1.fq',
        idx = expand(os.path.join(
            DIR_GRC_IDX, 'chr{}'.format(CHROM) + '_grc.{i}.bt2'), i = IDX_ITEMS)
    params:
        DIR_GRC_IDX + 'chr{}_grc'.format(CHROM)
    output:
        sam = os.path.join(DIR_FIRST_PASS, 'chr{}-grc.sam'.format(CHROM))
    threads: THREADS
    shell:
        '{BOWTIE2} --threads {threads} -x {params} -U {input.reads1} -S {output.sam}'

# rule align_to_per_haploid_setting:
#     input:
#         readsA1 = PREFIX_PER + '_hapA_1.fq',
#         readsB1 = PREFIX_PER + '_hapB_1.fq',
#         idxA = expand(DIR_PER_IDX + CHROM + '-per_hapA.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV),
#         idxB = expand(DIR_PER_IDX + CHROM + '-per_hapB.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
#     params:
#         indexA = DIR_PER_IDX + CHROM + '-per_hapA',
#         indexB = DIR_PER_IDX + CHROM + '-per_hapB'
#     output:
#         samA = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapA_haploid.sam'),
#         samB = os.path.join(DIR_FIRST_PASS, CHROM + '-per_hapB_haploid.sam')
#     threads: THREADS
#     shell:
#         '{BOWTIE2} --reorder --threads {threads} -x {params.indexA} -U {input.readsA1} -S {output.samA};'
#         '{BOWTIE2} --reorder --threads {threads} -x {params.indexB} -U {input.readsB1} -S {output.samB}'
# 
# PER_SAM_A = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapA.sam'.format(CHROM))
# PER_SAM_B = os.path.join(DIR_FIRST_PASS, 'chr{}-per_hapB.sam'.format(CHROM))
# 
# rule align_to_per:
#     input:
#         reads1 = PREFIX_PER + '_1.fq',
#         idxA = expand(DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapA.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV),
#         idxB = expand(DIR_PER_IDX + 'chr{}'.format(CHROM) + '-per_hapB.{IDX_ITEMS}.bt2', IDX_ITEMS = IDX_ITEMS, INDIV = INDIV)
#     params:
#         indexA = DIR_PER_IDX + 'chr{}-per_hapA'.format(CHROM),
#         indexB = DIR_PER_IDX + 'chr{}-per_hapB'.format(CHROM)
#     output:
#         samA = PER_SAM_A,
#         samB = PER_SAM_B
#     threads: THREADS
#     shell:
#         '{BOWTIE2} --reorder --threads {threads} -x {params.indexA} -U {input.reads1} -S {output.samA};'
#         '{BOWTIE2} --reorder --threads {threads} -x {params.indexB} -U {input.reads1} -S {output.samB}'
# 
# rule merge_per:
#     input:
#         samA = PER_SAM_A,
#         samB = PER_SAM_B
#     output:
#         path = os.path.join(DIR_FIRST_PASS, 'chr{}-per.paths'.format(CHROM)),
#         id = os.path.join(DIR_FIRST_PASS, 'chr{}-per.ids'.format(CHROM)),
#         merge_paths = os.path.join(DIR_FIRST_PASS, 'chr{}-per.merge_paths'.format(CHROM)),
#         samA = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA.sam'.format(CHROM)),
#         samB = os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB.sam'.format(CHROM))
#     params:
#         os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged'.format(CHROM))
#     run:
#         #: prepare ids
#         for h in ['hapA', 'hapB']:
#             shell('echo {h} >> {output.id};')
#         #: prepare paths
#         shell('ls {input.samA} >> {output.path};')
#         shell('ls {input.samB} >> {output.path};')
#         #: merge_incremental
#         shell('{PYTHON} {DIR_SCRIPTS}/merge_incremental.py -ns {output.path} \
#             -ids {output.id} -rs {RAND_SEED} -p {params} \
#             -l {output.merge_paths};')

rule check_standard_onepass:
    input:
        grc = expand(
            os.path.join(DIR_FIRST_PASS, 'chr{}-grc.sam'.format(CHROM)),
            INDIV = INDIV),
        # samA = expand(
        #     os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapA.sam'.format(CHROM)),
        #     INDIV = INDIV),
        # samB = expand(
        #     os.path.join(DIR_FIRST_PASS, 'chr{}-per-merged-hapB.sam'.format(CHROM)),
        #     INDIV = INDIV),
    output:
        touch(temp(os.path.join(DIR, 'standard_onepass.done')))
