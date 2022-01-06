# Summarize alignment results for real datasets
set -xp

MAPQ=10
LOG="30x.log"

print_summary () {
    printf "$1" >> $3
    printf "\t" >> $3
    printf `samtools view -c -F 256 -F 2048 $2` >> $3
    printf "\t" >> $3
    printf `samtools view -c -F 4 -F 256 -F 2048 $2` >> $3
    printf "\t" >> $3
    printf `samtools view -c -q ${MAPQ} -F 256 -F 2048 $2` >> $3
    printf "\n" >> $3
}

# Write header
printf "Sample\tAligner\tMethod\tNumReads\tNumMapped\tNumHighMapQ\n" > ${LOG}

# print_summary "HG002_0.3\tBWA-MEM\tCHM13" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg002_30x/0.3x_testdata/HG002-0.3x-bwa-chm13_v1.1_hg2y.bam ${LOG}
print_summary "HG002\tBWA-MEM\tGRCh38" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg002_30x/HG002.novaseq.pcr-free.30x.bwa.grch38.bam ${LOG}
print_summary "HG001\tBWA-MEM\tGRCh38" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg001_30x/HG001.novaseq.pcr-free.30x.bwa.grch38.bam ${LOG}
print_summary "HG005\tBWA-MEM\tGRCh38" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg005_30x/HG005.novaseq.pcr-free.30x.bwa.grch38.bam ${LOG}

print_summary "HG002\tBWA-MEM\tCHM13-to-GRCh38" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg002_30x/v0.5.1_new_lowmap/chm13_h38_exc_5/HG002-bwa-chm13_v1.1_hg2y_to_grch38-final.bam ${LOG}
print_summary "HG001\tBWA-MEM\tCHM13-to-GRCh38" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg001_30x/v0.5.1_new_lowmap/chm13_h38_exc_5/HG001-bwa-chm13_v1.1_hg2y_to_grch38-lowmap-unique_exc_1_5-final.bam ${LOG}
print_summary "HG005\tBWA-MEM\tCHM13-to-GRCh38" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg005_30x/v0.5.1_new_lowmap/chm13_h38_exc_5/HG005-bwa-chm13_v1.1_hg2y_to_grch38-lowmap-unique_exc_1_5-final.bam ${LOG}

print_summary "HG002\tBWA-MEM\tGRCh37" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg002_30x/HG002.novaseq.pcr-free.30x.bwa.grch37.bam ${LOG}
print_summary "HG001\tBWA-MEM\tGRCh37" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg001_30x/HG001.novaseq.pcr-free.30x.bwa.hg19.bam ${LOG}
print_summary "HG005\tBWA-MEM\tGRCh37" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg005_30x/HG005.novaseq.pcr-free.30x.bwa.hg19.bam ${LOG}

print_summary "HG002\tBWA-MEM\tCHM13-to-GRCh37" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg002_30x/v0.5.1_new_lowmap/chm13_h37_exc_1_5/HG002-bwa-chm13_v1.1_hg2y_to_hg19-lowmap-unique_exc_1_5-final.bam ${LOG}
print_summary "HG001\tBWA-MEM\tCHM13-to-GRCh37" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg001_30x/v0.5.1_new_lowmap/chm13_h37_exc_5/HG001-bwa-chm13_v1.1_hg2y_to_hg19-final.bam ${LOG}
print_summary "HG005\tBWA-MEM\tCHM13-to-GRCh37" /home/cnaechy1/scr16_blangme2/naechyun/leviosam_exp/hg005_30x/v0.5.1_new_lowmap/bwa_chm13_h37_exc_1_5/HG005-bwa-chm13_v1.1_hg2y_to_hg19-lowmap-unique_exc_1_5-final.bam ${LOG}
