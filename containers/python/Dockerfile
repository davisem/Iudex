# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Eric Davis" \
      description="Python Docker image" \
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

# Install Bedtools2, required for pybedtools.
RUN wget https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS2_VERSION/$BEDTOOLS2_SOURCE -O /opt/$BEDTOOLS2_SOURCE \
    && tar -xvf /opt/$BEDTOOLS2_SOURCE -C /opt \
    && cd /opt/bedtools2 \
    && make \
    && cp $BEDTOOLS2_BIN $BEDTOOLS2_DEST \
    && rm /opt/$BEDTOOLS2_SOURCE

# Finally, install python packages
RUN pip install -r ./requirements.txt
pip freeze > requirements.txt