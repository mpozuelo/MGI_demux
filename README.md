[![Build Status](https://travis-ci.org/mpozuelo/MGI_demux.svg?branch=master)](https://travis-ci.org/mpozuelo/MGI_demux)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/rnaseq.svg)](https://hub.docker.com/r/nfcore/rnaseq/)

### Introduction

**mpozuelo/MGI_demux** is a bioinformatics analysis pipeline used for demultiplexing Illumina data with i5 and i7 of.

The workflow processes raw data from FastQ inputs
([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),
 [cutadapt](https://github.com/marcelm/cutadapt/)) for demultiplexing
  See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Start running your own analysis!

You need to organize your fastq files as follows 'working_directory/data/fastq/run_id/lane/run_id_lane_read{1,2}.fq.gz'. Then from the working directory run the following instruction.

```bash
nextflow run mpozuelo/MGI_demux -profile <docker/singularity/conda> --input 'samplesheet.txt'
```

The samplesheet used must have at least 6 columns WITHOUT header: sample_ID, index_sequence, index2_sequence, barcode_sequence, run_ID, lane

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

### Documentation

The mpozuelo/MGI_demux pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Running the pipeline](docs/usage.md)

### Credits

These scripts were originally written for demultiplexing Illumina data with i5 and i7 adapters by Marta Pozuelo ([@mpozuelo](https://github.com/mpozuelo))
