# Performance

This pipeline has been assembled using the reactive workflow framework [Nextflow](https://www.nextflow.io/). Each process in the pipeline executed in a Docker container, because we're doing high impact science and reproducibility is important. Depending on your hardware, these design decisions could impact performance.

# OSX

There is currently no kernal support for Docker on OSX. This means Docker runs in a seperate kernal using virtualization and as a result, performance degrades severely, -especially when analyizing very large FASTQ files > 20G. Luckily you still have options for running *Iudex* at its peak performance.

**Option 1:**
Buy a Linux box.

**Option 2:**
Run `fastq_filterer` on bare metal. This will eleviate a significant ammount of the computation time in deduplicating your fastq files. To do this, compile `fastq_filterer.cpp` and run it outside of the pipeline to generate deduplicated fastqs. Then use these as inputs to the pipeline.

```
$ g++ -O3 fastq_filterer.cpp -o fastq_filterer
$ ./fastq_filterer duplicated.fastq output_name.fastq 0.0001
```


# Linux

Things are running optimally on Linux. 
