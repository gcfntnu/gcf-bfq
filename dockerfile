FROM continuumio/miniconda3

RUN conda update conda -y
RUN pip install --upgrade pip

RUN apt update && apt install -y -q git alien tabix pigz gzip openjdk-8-jdk openjdk-8-jre zlibc zlib1g-dev zlib1g
RUN conda update conda
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install -y fastq-screen bowtie2 fastqc bbmap pandas bioblend reportlab xlrd
RUN conda update numpy

WORKDIR /opt

RUN mkdir fastq_screen
COPY files/fastq_screen.conf fastq_screen/

RUN mkdir images
COPY files/ntnu.png images/

COPY files/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm /tmp
RUN alien -i /tmp/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm

RUN mkdir /bfq
WORKDIR /bfq
RUN mkdir output seq_share data log

WORKDIR /opt

RUN apt update && apt install -y make g++ libboost-iostreams-dev libboost-program-options-dev libz-dev
RUN git clone https://github.com/nsc-norway/suprDUPr
RUN cd suprDUPr && make 
ENV PATH=$PATH:/opt/suprDUPr

COPY files/cellranger-3.1.0.tar /tmp/cellranger-3.1.0.tar
RUN tar -C /opt/ -xvf /tmp/cellranger-3.1.0.tar
ENV PATH=$PATH:/opt/cellranger-3.1.0

COPY files/cellranger-atac-1.2.0.tar.gz /tmp/cellranger-atac-1.2.0.tar.gz
RUN tar -C /opt/ -xvf /tmp/cellranger-atac-1.2.0.tar.gz
ENV PATH=$PATH:/opt/cellranger-atac-1.2.0

COPY files/spaceranger-1.0.0.tar.gz /tmp/spaceranger-1.0.0.tar.gz
RUN tar -C /opt/ -xvf /tmp/spaceranger-1.0.0.tar.gz
ENV PATH=$PATH:/opt/spaceranger-1.0.0

RUN conda install -y configparser matplotlib
RUN conda install -y requests illumina-interop
RUN pip install cycler

RUN apt-get -qq update && apt-get install -y -q xkcdpass vim p7zip-full parallel
RUN conda install -y snakemake

WORKDIR /tmp

ENV DEBIAN_FRONTEND="noninteractive"
RUN wget -O- http://neuro.debian.net/lists/stretch.de-md.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN ["/bin/bash", "-c", "set -o pipefail && apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9"]
RUN apt-get update --yes && apt install -y singularity-container 

RUN ["/bin/bash", "-c", "set -o pipefail && curl -fsSL get.nextflow.io | /bin/bash"]
RUN mv nextflow /usr/local/bin

COPY files/logo_ntnu.png /opt/images/

ADD https://api.github.com/repos/gcfntnu/configmaker/git/refs/heads/master conf_version.json
RUN pip install --upgrade 'git+https://github.com/gcfntnu/configmaker#egg=configmaker'
ADD https://api.github.com/repos/flatberg/MultiQC/git/refs/heads/bfq mqc_version.json
#RUN pip install --upgrade 'git+https://github.com/flatberg/MultiQC/tree/bfq#egg=multiqc'
RUN pip install https://github.com/flatberg/MultiQC/archive/bfq.zip

#FOR DEV
#ADD https://api.github.com/repos/gcfntnu/gcfdb/git/refs/heads/dev conf_version.json
#RUN cd /opt && git clone https://github.com/gcfntnu/gcfdb.git && cd /opt/gcfdb && git checkout dev 
RUN cd /opt && git clone https://github.com/gcfntnu/gcfdb.git && cd /opt/gcfdb && git checkout tags/0.1
ENV GCF_DB=/opt/gcfdb

#FOR DEV
#ADD https://api.github.com/repos/gcfntnu/rna-seq/git/refs/heads/dev conf_version.json
#RUN cd /opt && git clone https://github.com/gcfntnu/rna-seq.git && cd /opt/rna-seq && git checkout dev
RUN cd /opt && git clone https://github.com/gcfntnu/rna-seq.git && cd /opt/rna-seq && git checkout tags/0.1

COPY ./bcl2fastq_pipeline /opt/bcl2fastq_pipeline
ENV PATH=$PATH:/opt/bcl2fastq_pipeline/flowcell_manager

WORKDIR /opt/bcl2fastq_pipeline
RUN gcc splitFastq.c -o splitFastq
RUN ./setup.py build
RUN ./setup.py install

ENV PATH=$PATH:/opt/bcl2fastq_pipeline

ADD ./files/snakefiles /opt/snakefiles

ENV TMPDIR=/instrument-archive/tmp
ENV EXT_DIR=/instrument-archive/ext
ENV GCF_EXT=$EXT_DIR
ENV SINGULARITY_BINDPATH=/mnt,/instrument-archive,/opt/gcfdb
ENV SINGULARITY_CACHEDIR=$TMPDIR/singularity/
ENV SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
ENV SINGULARITY_LOCALCACHEDIR=$SINGULARITY_CACHEDIR/tmp/
ENV SINGULARITY_PULLFOLDER=$SINGULARITY_CACHEDIR

ENV PATH=$PATH:/opt/bcl2fastq_pipeline
CMD bin/bfq.py
