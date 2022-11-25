from collections import defaultdict
import pysam
import sys

with pysam.VariantFile(sys.argv[1]) as f:
    chr_pos = set()
    cat_cnt = defaultdict(int)
    for record in f:
        if (record.chrom, record.start, record.alleles) in chr_pos:
            continue
        chr_pos.add((record.chrom, record.start, record.alleles))
        cat_cnt[record.info['SVTYPE']] += 1

print(cat_cnt)