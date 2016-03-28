#/Groups/twntyfr/bin/bowtie2-2.1.0/bowtie2 --sensitive -x /Groups/twntyfr/index/hg19 -U 001A_ATCACG_L001_R1_001.fastq -S test.sam --un unaligned.out -p 2 & 



SAMPLE="020B"

nohup /Groups/twntyfr/bin/bowtie2-2.1.0/bowtie2 --very-sensitive -x /Groups/twntyfr/index/hg19 -U $SAMPLE.all.fastq -S $SAMPLE.human.sam --un $SAMPLE.nonhuman.fastq -p 8 2> $SAMPLE.human.summary.txt & 




