#!/bin/bash
set -x #echo on

#
# Example: kallisto.sh /u/eaudemard/project/epcy_paper/data/ SRR2753073
#

disease=$1
sample=$2

module load kallisto/0.43.1
module load sratoolkit

DESTDIR=$1/$sample
WORKDIR=$TMPDIR/$sample

mkdir -p $WORKDIR
mkdir -p $DESTDIR

fastq-dump --split-files $sample -O $WORKDIR && \
    java -Xms4g -Xmx4g -jar /soft/bioinfo/linux_RH7/mugqic_pipelines-2.2.0/resources/software/trimmomatic/Trimmomatic-0.35/trimmomatic-0.35.jar PE \
    		     -threads 8 -phred33 \
    		     $WORKDIR/${sample}_1.fastq \
    		     $WORKDIR/${sample}_2.fastq \
    		     $WORKDIR/${sample}_1.ffastq \
    		     $WORKDIR/${sample}_1.unpaired.fq \
    		     $WORKDIR/${sample}_2.ffastq \
    		     $WORKDIR/${sample}_2.unpaired.fq \
		     ILLUMINACLIP:/soft/bioinfo/linux_RH7/mugqic_pipelines-2.2.0/resources/software/trimmomatic/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10:8:true \
    		     LEADING:20 TRAILING:20 MINLEN:20 && \
    ls $WORKDIR && \
    kallisto quant -b 20 -t 8 -i /share/tsa_project/lib/kallisto/GRCh38.88/GRCh38.88.cdna.kai -o $DESTDIR/kallisto.GRCh38.Gencode88.transcripts $WORKDIR/${sample}_1.ffastq $WORKDIR/${sample}_2.ffastq && \
    touch $DESTDIR/$sample.GRCh38.Gencode88.transcripts.done
