FROM continuumio/miniconda3

COPY environment.yml .

RUN conda env create -f environment.yml

SHELL ["conda", "run", "-n", "fastq-pipeline", "/bin/bash", "-c"]

ENV PATH=/opt/conda/envs/fastq-pipeline/bin:$PATH
