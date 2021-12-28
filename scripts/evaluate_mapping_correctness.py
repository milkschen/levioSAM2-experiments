import pysam
import sys

TOLERANCE = 10

gold = {}
with pysam.AlignmentFile(sys.argv[1], 'r') as f:
    for r in f:
        if r.is_read1:
            gold[r.query_name + '_1'] = (r.reference_name, r.reference_start)
        elif r.is_read2:
            gold[r.query_name + '_2'] = (r.reference_name, r.reference_start)
print(len(gold))

cnt_correct = 0
cnt_all = 0
with pysam.AlignmentFile(sys.argv[2], 'r') as f:
    for r in f:
        if r.is_secondary or r.is_supplementary:
            continue

        cnt_all += 1
        if r.is_read1:
            g = gold[r.query_name + '_1']
            if g[0] == r.reference_name and abs(g[1] - r.reference_start) <= TOLERANCE:
                cnt_correct += 1
        elif r.is_read2:
            g = gold[r.query_name + '_2']
            if g[0] == r.reference_name and abs(g[1] - r.reference_start) <= TOLERANCE:
                cnt_correct += 1

print(cnt_all)
print(cnt_correct)
print(cnt_correct / cnt_all)
