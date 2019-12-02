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

RUN conda install -y configparser matplotlib
RUN conda install -y requests illumina-interop
RUN pip install cycler

RUN apt-get -qq update && apt-get install -y -q xkcdpass vim p7zip-full parallel

COPY files/logo_ntnu.png /opt/images/

ADD https://api.github.com/repos/gcfntnu/configmaker/git/refs/heads/master version.json
RUN pip install --upgrade 'git+https://github.com/gcfntnu/configmaker#egg=configmaker'
ADD https://api.github.com/repos/flatberg/MultiQC/git/refs/heads/master version.json
RUN pip install --upgrade 'git+https://github.com/flatberg/MultiQC#egg=multiqc'

#RUN conda install -y xlrd>=1.0.0
#ADD https://api.github.com/repos/gsmashd/bcl2fastq_pipeline/git/refs/heads/master version.json
#ARG CACHE_BREAK=random_num
#RUN git clone https://github.com/gsmashd/bcl2fastq_pipeline.git
COPY ./bcl2fastq_pipeline /opt/bcl2fastq_pipeline
ENV PATH=$PATH:/opt/bcl2fastq_pipeline/flowcell_manager

WORKDIR /opt/bcl2fastq_pipeline
RUN pwd && ls -la
RUN gcc splitFastq.c -o splitFastq
RUN ./setup.py build
RUN ./setup.py install

ENV PATH=$PATH:/opt/bcl2fastq_pipeline

CMD bin/bfq.py
