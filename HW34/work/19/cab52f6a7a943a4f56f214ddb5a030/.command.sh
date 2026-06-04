#!/bin/bash -ue
bwa index yeast_genome.fasta

bwa mem yeast_genome.fasta SRR6357070_1_trimmed.fq.gz             | samtools sort -o aligned.bam

samtools index aligned.bam
