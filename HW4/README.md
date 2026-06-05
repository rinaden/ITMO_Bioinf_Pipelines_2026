To run pipeline use this command:

```
nextflow run main.nf \
  -profile local \
  --samplesheet samples.csv \
  --reference ecoli_reference.fna
```

Stub mode run:

```
nextflow run main.nf \
  -profile local \
  -stub-run \
  --samplesheet samples.csv \
  --reference ecoli_reference.fna
```

Profile could be either local, container or cluster
