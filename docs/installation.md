Installation
------------

### Step 1 -- Download Nextflow
```
$ curl -fsSL get.nextflow.io | bash
$ ./nextflow
$ mv nextflow /usr/local/bin
```

### Step 2 -- Download Docker Containers

This `docker pull` command will download each of the required tools to run the pipeline.
```
$ docker pull emdavis/iudex -a 
```
The following images should be pulled down from Docker Hub.
```
bwa: Pulling from emdavis/iudex
Digest: sha256:8f5d1d70831f8a52fb92bc67a10a3bd5928f11953d589f386e0cca97e0fbb76e
fastq_filterer: Pulling from emdavis/iudex
Digest: sha256:3c68785e7e2dbd57f8acc9a6d481ed03a8b60d71a7dfe1f65b0060b78b4ccccd
insertion_counts: Pulling from emdavis/iudex
Digest: sha256:e65aa5148f0555e0d019eccf9acfa1d3dbce9cf0818407cf4463f28f24c0467c
Status: Image is up to date for emdavis/iudex
```

### Step 3 -- Download Github Repository
```
$ git clone https://github.com/davisem/Iudex
$ cd Iudex/
```

### Step 4 -- Create a BWA index
```
$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O annotations/hg38
$ cd annotations/hg38
$ gunzip hg38.fa.gz
$ docker run edavis/bwa:1.0 /bin/sh bwa index hg38.fa
```
### Step 5 -- Download Iudex repo
```
$ git clone https://github.com/davisem/Iudex.git
$ cd Iudex
```

### Step 6 -- Run some test data
```
$ nextflow iudex.nf --output_dir "/your/output/directory" --index "path/to/your/genome/index/folder/*"
```
