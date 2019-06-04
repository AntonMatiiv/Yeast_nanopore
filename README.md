# Analysis of yeast genomes from the Peterhof genetic collection
## Description
The Peterhof genetic collection of Saccharomyces cerevisiae strains (PGC) is a large laboratory fund, which has accumulated several thousand strains for more than half a century. Several PGC strains have been widely used in certain areas of yeast research, but their genomes have not yet been fully studied. The genetic distance between the precursor PGC and S288C is comparable to that between two geographically isolated populations.This project is a continuation of a project to assemble the yeast genome from Oxford Nanopore data. During this project it is supposed to complete the assembly of the reference genome of strain 1A-D1628 from PGC, and also to compare the genomes of other strains (74-D694, etc.) with the obtained reference assembly. Also, it is important to find genetic variants associated with certain features of the phenotype of mutant derivatives 1A-D1628 (based on data obtained using the Illumina technology).

## Goals and objectives
### Goal

The goal of the project is to identify compensatory mutations that allow cells to survive with the disruption of the gene encoding the vital translation termination factor.

### Objectives

1. To polish the assembly of genome 1A-D1628 with raw Oxford Nanopore technologies (ONT) data and Illumina readings.
2. To carry out a comparative analysis of the resulting assembly with S288C and strain 74-D694, sequenced on the Oxford Nanopore technologies platform.
3. To analyze the spectrum of genetic variants in the genomes of mutant derivatives 1A-D1628, capable of maintaining viability against the background of translation termination defects (using data from the sequencing of the mutant genomes on the Illumina platform).

## Requirements

* nanopolish version 0.11.1
* samtools version 1.9
* minimap2 version 2.14-r894
* QUAST version 4.5
* racon version v1.3.3
* exonerate version 2.4.0
* snpEff version 4.3t
* GATK version 4.1.2.0
* Picard version 2.0.1
* seqkit version v0.10.1
* bwa version 0.7.17-r1188

## Files description

* ```report.pdf```  - report output file of QUAST
* ```calling.sh``` - script for variant calling with the GATK
* ```QUAST_pics``` - images from QUAST report
* ```AC_F.jpg``` - allele counts frequency histogram

## Main commands
### Duplicates removal

Firstly, we need to remove duplicates from ONT reads by using ``seqkit`` with ONT reads:

```seqkit rmdup yeast_1d.guppy_213.fastq.gz > yeast_rmdup.fastq.gz```

Where ```yeast_1d.guppy_213.fastq.gz``` is ONT reads and ```yeast_rmdup.fastq.gz``` is ONT reads without duplicates.

### Genome polishing
#### Data preprocessing

The original purpose of nanopolish is to improve the consensus accuracy of an assembly of Oxford Nanopore Technology sequencing reads. Nanopolish needs access to the signal-level data measured by the nanopore sequencer. To begin, we need to create an index readdb file that links read ids with their signal-level data in the FAST5 files:

```nanopolish index -d ./uncleaned_fast5 yeast_rmdup.fastq.gz```

Where ```./uncleaned_fast5``` is the path to FAST5 files.

We get the following files: ```yeast_rmdup.fasta.gz.index```, ```yeast_rmdup.fasta.gz.index.fai```, ```yeast_rmdup.fasta.gz.index.gzi```, and ```yeast_rmdup.fasta.gz.index.readdb```.

#### Compute a new consensus sequence for a draft assembly

Now that we have ```yeast_rmdup.fastq.gz``` indexed with ```nanopolish index```, and have a previously assembled draft genome ```yeast_1d.guppy_213.canu.fa```, we can begin to improve the assembly with nanopolish.

First, we align the original reads (```yeast_rmdup.fastq.gz```) to the draft assembly (```yeast_1d.guppy_213.canu.fa```) and sort alignments:

```minimap2 -ax map-ont -t 8 yeast_1d.guppy_213.canu.fa yeast_rmdup.fastq.gz | samtools sort -o read_after_minimap2_sorted.bam -T reads.tmp``` 
```samtools index read_after_minimap2_sorted.bam```

Then we run the consensus algorithm. We use ```nanopolish_makerange.py``` to split the draft genome assembly into 50kb segments, so that we can run the consensus algorithm on each segment in parallel:

```python nanopolish_makerange.py yeast_1d.guppy_213.canu.fa | parallel --results nanopolish.results -P 2 nanopolish/nanopolish variants --consensus -o polished.{1}.vcf -w {1} -r yeast_rmdup.fasta.gz -b read_after_minimap2_sorted.bam -g yeast_1d.guppy_213.canu.fa -t 4 --min-candidate-frequency 0.1 ```

After all polishing jobs are complete, me can merge the individual 50kb segments together back into the final assembly:

```nanopolish vcf2fasta -g yeast_1d.guppy_213.canu.fa polished.*.vcf > polished_genome.fa```

### Racon polishing

Racon can be used as a polishing tool after the assembly with either Illumina data or data produced by third generation of sequencing. The type of data inputed is automatically detected. Racon takes as input only three files: contigs in FASTA/FASTQ format, reads in FASTA/FASTQ format and overlaps/alignments between the reads and the contigs in MHAP/PAF/SAM format. 

```bwa index  polished_genome.fa```

Where ``` polished_genome.fa``` is polished genome with ONT reads.

```bwa mem polished_genome.fa M1628_PE_merged.fastq.gz > alignment_polished_wt_PE_MERGED.sam```

Where ```M1628_PE_merged.fastq.gz``` are Illumina reads.

```racon M1628_PE_merged.fastq.gz alignment_polished_wt_PE_MERGED.sam polished_genome.fa > polished_genome_racon.fa```

### QUAST evaluation of genome assemblies

QUAST evaluates genome assemblies by computing various metrics. It works both with and without reference genomes. The tool accepts multiple assemblies, thus is suitable for comparison.

```quast -R S288C.ref.fa yeast_1d.guppy_213.canu.fa polished_genome.fa polished_genome_racon.fa```

Where ```S288C.ref.fa``` is S288C strain reference genome.

As a result we obtain some metrics calculated:



### Functional annotation of polished genome

Exonerate is a generic tool for pairwise sequence comparison. It allows you to align sequences using a many alignment models. We use ```protein2genome``` model and we will report the best 1 result for each protein query:

```exonerate --model protein2genome --bestn 1 uniprot-proteome%3AUP000002311.fasta polished_genome_racon.fa --showtargetgff yes --showvulgar no --showalignment no > exonerate_p2g_bestn1.gtf```

Where ```uniprot-proteome%3AUP000002311.fasta``` list of protein amino acid sequences of S288C from UniProt database (https://www.uniprot.org/proteomes/UP000002311) and ```exonerate_p2g_bestn1.gtf``` annotation in gff2 format.

### Building Snpeff Database

SnpEff needs a database to perform genomic annotations. In order to build a database for a new genome, we need to: 

```
export DBNAME=Saccharomyces_cerevisiae_1A-D1628_wt
export GFF=/path/to/exonerate_p2g_bestn1.gtf
export FASTA =/path/to/polished_genome_racon.fa

#Go into the snpEff directory and create a directory for your files
cd /usr/local/snpEff
mkdir data/$DBNAME

#Copy the files into snpEff's directory structure
cp $GFF data/$DBNAME/genes.gff
cp $FASTA data/$DBNAME/sequences.fa
 
#Edit snpEff.config and insert your specific database information:
echo "$DBNAME.genome : $DBNAME" >> snpEff.config

#Build the database
java -jar snpEff.jar build -gff2 -v $DBNAME
```
 
 ### Variant calling with the GATK
 
 Calling SNPs and indels via local re-assembly of haplotypes with ```calling.sh``` script.
 
 Then perform joint genotyping on gVCF files produced by HaplotypeCaller:
 
 ```
samples=$(find . | sed 's/.\///' | grep -E 'g.vcf$' | sed 's/^/--variant /')
java -Xmx8g -jar /path/to/GATK -T GenotypeGVCFs -R /path/to/polished_genome_racon.fa -o $(echo $samples) > M_total_ploidy1.g.vcf
```
 
 ### Examination of allele counts in g.vcf for mutants
 
Bash:
```
cut -f 8 M_total_ploidy1.g.vcf | cut -c 1-6 | grep -o '^AC=[0-9]*' > ac
grep -o '[0-9]*' ac > numbers_ac
```
Where ```M_total_ploidy1.g.vcf``` is the result of ```GenotypeGVCFs```

R:
```
numbersac <- read.table("/path/to/numbers_ac", header = TRUE)
hist(numbersac$X1, freq = TRUE, col = "blue", xlab="Allele counts", ylab="Frequency", main = "", nclass = 100)
```


