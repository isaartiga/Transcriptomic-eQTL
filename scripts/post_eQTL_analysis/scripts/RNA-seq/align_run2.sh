#!/bin/bash

PROJECT_HOME=/mnt/sdb/projects/Iraide/MS_Iraide_Alloza_2023_10/RNA_seq/fastq/
PACKAGE=$PROJECT_HOME/20230616_Iraide_RNAseq
REFERENCE_INDEX=/mnt/sdb/references/Genomes/Hg38/Gencode
OUTPUT=/mnt/sdb/projects/Iraide/MS/RNA_seq/bams


CPUS=8

#. /etc/profile.d/lmod.sh
module load subread/2.0.6

for i in "$PACKAGE/"*_R1_001.fastq.gz; do
  TRIMMED_P1="$i"
  TRIMMED_P2="$(echo "$i" | sed 's/_R1/_R2/g')"
  STAR_MAPPING_PREFIX=`basename $TRIMMED_P1 | sed s/_R1_001\.fastq\.gz//`

  echo "$TRIMMED_P1"
  echo "$TRIMMED_P2"
  echo "$STAR_MAPPING_PREFIX"
subread-align -t 0 -T 60  -i $REFERENCE_INDEX/GRCh38 -r $TRIMMED_P1 -R $TRIMMED_P2 -o $OUTPUT/$STAR_MAPPING_PREFIX.bam

done
