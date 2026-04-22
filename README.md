# ITMO_Bioinf_Pipelines_2026
Repository for homeworks

**Comments for HW2**
To run pipeline you can use this command:

```
nextflow run main.nf \
  --reads "test_data/SRR6357070_{1,2}.fastq.gz" \
  --reference test_data/yeast_genome.fasta \
  --paired true \
  --outdir results \
  -with-conda
```

Pipeline accepts both single- and paired-end reads, depending on the `--paired` flag. Use "true" for paired-end and "false" for single-end. 
