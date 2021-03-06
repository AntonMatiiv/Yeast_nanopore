#!/bin/bash

export GATK='/path/to/gatk-protected/target/executable/GenomeAnalysisTK.jar'
export PICARD='/path/to/picard-tools-2.0.1/picard.jar'
mkdir bams gvcfs logs
for i in *_R1_*fastq.gz
#Where *_R1_*fastq.gz is regex for mutants Illumina reads
do
        bwa mem -t 32 -R "@RG\tID:${i%%_R1*}\tSM:S${i%%_*}\tLB:1\tPL:illumina" $i ${i%%_R1*}_R2*gz | /path/to/software/samtools view -bS - > /path/to/bams/${i%%_R1*}.bam
done
cd bams
for i in M*bam
#Where M*bam is bam files for mutants
do
        /path/to/samtools sort -T ${i%%.bam} -o ${i%%.bam}.sorted.bam $i &
done
wait
for i in M*.sorted.bam
#Where M*.sorted.bam is regex for mutants
is bam files for mutants
do
        mv $i ${i%%_L008*}.bam
done
for i in *L001*bam
do
       while [ $( ps -u $USER | grep 'samtools' | wc -l ) -ge 24 ] ; do sleep 1 ; done
       /path/to/samtools merge ${i%%_L001*}.bam ${i%%L001*}*bam &
done
wait
rm *L00*bam
mkdir tmp
for BAM in M*.bam
do
        java -Xmx2g -jar $PICARD MarkDuplicates I=$BAM O=${BAM%%.bam}.dedup.bam M=/path/to/logs/${BAM%%.bam}.MD.metrics ASSUME_SORTED=true TMP_DIR=${PWD}/tmp &
done
wait
rm -rf tmp
for i in M*.dedup.bam
#Where M*.dedup.bam is regex for mutants
do
       samtools index $i &
done
wait
for i in M*.dedup.bam
do
        java -Xmx8g -jar $GATK -T HaplotypeCaller -R /path/to/polished_genome_racon.fa -ERC GVCF -ploidy 1 -I $i -o /path/to/gvcfs/${i%%.ded*}.raw.g.vcf &
done

