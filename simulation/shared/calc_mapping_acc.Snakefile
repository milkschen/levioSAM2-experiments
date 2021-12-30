'''
Calculate mapping accuracy for all settings
'''
rule calc_mapping_accuracy_target:
    input:
        truth = os.path.join(PREFIX_PER + '.bam'),
        exp = os.path.join(
            DIR_FIRST_PASS,
            '{aln}' + f'-{TARGET_LABEL}.bam'),
    output:
        acc = os.path.join(
            DIR_RESULTS,
            '{INDIV}-{aln}' + f'-{TARGET_LABEL}.acc')
    shell:
        '{PYTHON} {DIR_SCRIPTS}/evaluate_mapping_correctness.py -g {input.truth} -q {input.exp} '
        '-o {output.acc}'

rule calc_mapping_accuracy_leviosam:
    input:
        truth = os.path.join(PREFIX_PER + '.bam'),
        exp = os.path.join(
            DIR_FIRST_PASS,
            '{aln}' + f'-{SOURCE_LABEL}_to_{TARGET_LABEL}-final.bam'),
    output:
        acc = os.path.join(
            DIR_RESULTS,
            '{INDIV}-{aln}' + f'-{SOURCE_LABEL}_to_{TARGET_LABEL}.acc')
    shell:
        '{PYTHON} {DIR_SCRIPTS}/evaluate_mapping_correctness.py -g {input.truth} -q {input.exp} '
        '-o {output.acc}'

'''
Summarize results as a TSV
'''
rule check_mapping_acc_and_write_as_tsv:
    input:
        acc = expand(os.path.join(
            DIR_RESULTS,
            '{indiv}-{aln}' + '-{method}.acc'),
            indiv = INDIV, aln = ['bt2', 'bwa'],
            method=[TARGET_LABEL, f'{SOURCE_LABEL}_to_{TARGET_LABEL}']
        )
    output:
        tsv = os.path.join(DIR_RESULTS, 'all.tsv'),
        check = touch(temp(os.path.join(DIR, 'accuracy.done')))
    run:
        dict_indiv_to_pop = build_dict_indiv_to_pop(FAMILY)
        dict_pop_to_spop = build_dict_pop_to_spop(SPOP)
        
        df = pd.DataFrame()
        list_indiv = []
        list_aln = []
        list_method = []
        list_all = []
        list_mapped = []
        list_high_mapq = []
        list_tp = []
        for fn in input.acc:
            if fn.endswith('.acc'):
                d = {}
                with open(fn, 'r') as f:
                    for line in f:
                        line = line.split()
                        d[line[0]] = line[1]
                    list_all.append(int(d['num_query']))
                    list_mapped.append(int(d['num_mapped']))
                    list_high_mapq.append(int(d['num_high_mapq']))
                    list_tp.append(int(d['num_correct']))
                fn = fn[: fn.rfind('.acc')]
                bn = os.path.basename(fn).split('-')
                list_indiv.append(bn[0])
                list_aln.append(bn[1])
                list_method.append(bn[2])
        df['Inidvidual'] = list_indiv
        df['SuperPopulation'] = [dict_pop_to_spop[dict_indiv_to_pop[s]] for s in list_indiv]
        df['Method'] = list_method
        df['Aligner'] = list_aln
        df['NumReads'] = list_all
        df['NumMapped'] = list_mapped
        df['NumHighMapq'] = list_high_mapq
        df['TruePositive'] = list_tp
        df.to_csv(output.tsv, sep='\t', index=None, float_format = '%.4f')
