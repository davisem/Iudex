# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Eric Davis" \
      description="BWA (v0.7.15) Docker image" \
      maintainer="emdavis48@gmail.com"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
 && rm -rf /var/lib/apt/lists/*

COPY fastq_filterer.cpp /src/
COPY bloom/bloom_filter.hpp /src/

RUN g++ -O2 /src/fastq_filterer.cpp -o /src/fastq_filterer
CMD [./fastq_filterer]