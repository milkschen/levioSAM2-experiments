'''
Simple mapping accuracy evaluation of two SAM/BAM files

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import pysam
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g', '--gold', required=True,
        help='Path to the truth SAM/BAM file. [required]'
    )
    parser.add_argument(
        '-q', '--query', required=True,
        help='Path to the query SAM/BAM file. [required]'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output report. Leave empty to write to stdout. [empty]'
    )
    parser.add_argument(
        '-t', '--tolerance', default=10,
        help='Allowed bases of position difference (diff > `-t` is consider incorrect). [10]'
    )
    parser.add_argument(
        '-m', '--high_mapq', default=10,
        help='Lower bound MAPQ for an alignment to be considered as high-quality. [10]'
    )
    args = parser.parse_args()
    return args


def evaluate_mapping_correctness(fn_gold, fn_query, f_out, tolerance, high_mapq):
    gold = {}
    with pysam.AlignmentFile(fn_gold, 'r') as f:
        for r in f:
            if r.is_read1:
                gold[r.query_name + '_1'] = (r.reference_name, r.reference_start)
            elif r.is_read2:
                gold[r.query_name + '_2'] = (r.reference_name, r.reference_start)

    cnt_correct = 0
    cnt_mapped = 0
    cnt_high_mapq = 0
    cnt_all = 0
    printed_err = 0
    with pysam.AlignmentFile(fn_query, 'r') as f:
        for r in f:
            if r.is_secondary or r.is_supplementary:
                continue

            cnt_all += 1
            if not r.is_unmapped:
                cnt_mapped += 1
            if r.mapping_quality >= high_mapq:
                cnt_high_mapq += 1

            if r.is_read1:
                g = gold[r.query_name + '_1']
                if g[0] == r.reference_name and abs(g[1] - r.reference_start) <= tolerance:
                    cnt_correct += 1
                elif printed_err <= 30:
                    # print(g)
                    # print(r)
                    printed_err += 1
            elif r.is_read2:
                g = gold[r.query_name + '_2']
                if g[0] == r.reference_name and abs(g[1] - r.reference_start) <= tolerance:
                    cnt_correct += 1

    print(f'num_gold\t{len(gold)}', file=f_out)
    print(f'num_query\t{cnt_all}', file=f_out)
    print(f'num_mapped\t{cnt_mapped}', file=f_out)
    print(f'num_high_mapq\t{cnt_high_mapq}', file=f_out)
    print(f'num_correct\t{cnt_correct}', file=f_out)
    print(f'correct/all\t{cnt_correct / cnt_all:.6f}', file=f_out)
    print(f'correct/mapped\t{cnt_correct / cnt_mapped:.6f}', file=f_out)

    return


if __name__ == '__main__':
    args = parse_args()
    fn_gold = args.gold
    fn_query = args.query
    tolerance = args.tolerance
    high_mapq = args.high_mapq

    if args.output != '':
        f_out = open(args.output, 'w')
    else:
        f_out = sys.stdout

    evaluate_mapping_correctness(
        fn_gold=fn_gold, fn_query=fn_query, f_out=f_out,
        tolerance=tolerance, high_mapq=high_mapq
    )
