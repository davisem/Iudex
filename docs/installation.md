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

### Step 3 -- Download Github Repository
```
$ git clone https://github.com/davisem/Iudex
$ cd Iudex/
```
