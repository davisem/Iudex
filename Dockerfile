# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Eric Davis" \
      description="Iudex Docker image" \
      maintainer="emdavis48@gmail.com"

# Install core libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    python \
    python-pip \
    wget \
    build-essentials \
    zlib \
 && rm -rf /var/lib/apt/lists/*


# Install Bedtools2
RUN wget https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS2_VERSION/$BEDTOOLS2_SOURCE -O /opt/$BEDTOOLS2_SOURCE \
    && tar -xvf /opt/$BEDTOOLS2_SOURCE -C /opt \
    && cd /opt/bedtools2 \
    && make \
    && cp $BEDTOOLS2_BIN $BEDTOOLS2_DEST \
    && rm /opt/$BEDTOOLS2_SOURCE

# Install Pandas
RUN pip install pybedtools \
    pip install pandas

# Install Bowtie2
RUN git clone https://github.com/BenLangmead/bowtie2.git && \
    cd bowtie2 && \
    make NO_TBB=1 && \
    make install && \
    cd && \
    rm -rf /bowtie2

#Include bed annotations
COPY annotations/*