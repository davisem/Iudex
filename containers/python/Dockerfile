# Set the base image
FROM debian:jessie-slim

# Set image metadata
LABEL author="Eric Davis" \
      description="Python Docker image" \
      maintainer="emdavis48@gmail.com"

RUN apt-get update -y && apt-get install
RUN apt-get install -y python-dev python-pip

RUN pip install -r ./requirements.txt
pip freeze > requirements.txt