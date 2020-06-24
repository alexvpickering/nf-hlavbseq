FROM openjdk:8-jre
LABEL maintainer="Alex Pickering <alexvpickering@gmail.com>"

ENV BWA_VERSION 0.7.17
ENV PICARD_VERSION 2.18.21
ENV SAMTOOLS_VERSION 1.9

WORKDIR /tmp

RUN apt-get update -y \
  && apt-get install --no-install-recommends -y \
       make \
       gcc \
       g++ \
       libz-dev \
       libbz2-dev \
       liblzma-dev \
       ncurses-dev \
       libnss-sss \
       time \
  && wget -q https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 \
  && tar xjvf bwa-${BWA_VERSION}.tar.bz2 \
  && cd /tmp/bwa-${BWA_VERSION}/ \
  && make \
  && cp -av /tmp/bwa-${BWA_VERSION}/bwa /usr/bin/ \
  && cd /tmp \
  && wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
  && cd /tmp/samtools-${SAMTOOLS_VERSION}/ \
  && make \
  && cp -av /tmp/samtools-${SAMTOOLS_VERSION}/samtools /usr/bin/ \
  && wget -q -O /usr/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
  && ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime \
  && echo "America/Chicago" > /etc/timezone \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get clean all \
  && rm -rfv /var/lib/apt/lists/* /tmp/* /var/tmp/* 

RUN mkdir ~/hg19 && cd ~/hg19 \
    && wget https://storage.googleapis.com/genomics-public-data/references/b37/Homo_sapiens_assembly19.fasta.gz

RUN cd ~/hg19 && gunzip Homo_sapiens_assembly19.fasta.gz \
    && bwa index Homo_sapiens_assembly19.fasta

RUN mkdir ~/HLAVBSeq && cd ~/HLAVBSeq \
    && wget http://nagasakilab.csml.org/hla/HLAVBSeq.jar \
    && wget http://nagasakilab.csml.org/hla/bamNameIndex.jar \
    && wget http://nagasakilab.csml.org/hla/parse_result.pl \
    && wget http://nagasakilab.csml.org/hla/call_hla_digits.py \
    && wget http://nagasakilab.csml.org/hla/hla_all_v2.fasta \
    && wget http://nagasakilab.csml.org/hla/Allelelist_v2.txt

COPY --from=python:2.7-alpine3.11 / /