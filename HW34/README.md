To run pipeline use this command:

```
nextflow run main.nf \
  -profile local \
  --reads test_data/SRR6357070_{1,2}.fastq.gz \
  --reference test_data/yeast_genome.fasta \
  --paired true \
  --outdir results
```

Profile could be either local, container or cluster
