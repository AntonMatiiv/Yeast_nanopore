# Analysis of yeast genomes from the Peterhof genetic collection
## Description
The Peterhof genetic collection of Saccharomyces cerevisiae strains (PGC) is a large laboratory fund, which has accumulated several thousand strains for more than half a century. Several PGC strains have been widely used in certain areas of yeast research, but their genomes have not yet been fully studied. The genetic distance between the precursor PGC and S288C is comparable to that between two geographically isolated populations.This project is a continuation of a project to assemble the yeast genome from Oxford Nanopore data. During this project it is supposed to complete the assembly of the reference genome of strain 1A-D1628 from PGC, and also to compare the genomes of other strains (74-D694, etc.) with the obtained reference assembly. Also, it is important to find genetic variants associated with certain features of the phenotype of mutant derivatives 1A-D1628 (based on data obtained using the Illumina technology).

## Goals
1. To polish the assembly of genome 1A-D1628 with raw Oxford Nanopore technologies (ONT) data and Illumina readings.
2. To carry out a comparative analysis of the resulting assembly with S288C and strain 74-D694, sequenced on the Oxford Nanopore technologies platform.
3. To analyze the spectrum of genetic variants in the genomes of mutant derivatives 1A-D1628, capable of maintaining viability against the background of translation termination defects (using data from the sequencing of the mutant genomes on the Illumina platform).

## Requirements

* nanopolish
* samtools
* minimap2
* QUAST
* racon
* exonerate
* snpEff
* GATK
* seqkit

## Files description

* ```main_commands``` - project main commands and their description
* ```report.pdf```  - report output file of QUAST 
* yeast_1d.guppy_213.fastq.gz - ONT reads
* 
