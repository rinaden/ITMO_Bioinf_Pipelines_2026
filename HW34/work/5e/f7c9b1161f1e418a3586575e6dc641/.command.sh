#!/bin/bash -ue
samtools index aligned.bam

bcftools mpileup         -f yeast_genome.fasta         aligned.bam     | bcftools call         -mv -Oz         -o SRR6357070_1.vcf.gz

tabix -p vcf SRR6357070_1.vcf.gz
