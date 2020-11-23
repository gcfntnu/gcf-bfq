#!/bin/bash

if [ "$#" -lt 2 ] ; 
then
    echo  -e "Usage: \n" \
    "./build-tag-push.sh [prod|test|base] [docker-tag] [optional args]\n\n" \
    "Optional branch args for dockerfile-test:\n" \
    "-d gcfdb \n" \
    "-r rna-seq \n" \
    "-s single-cell \n" \
    "-t gcf-tools \n" \
    "-q multiqc \n" \
    "-m microbiome \n" 
	exit 1
fi

gcf_db=bfq-dev
rna_seq=bfq-dev
single_cell=bfq-dev
gcf_tools=bfq-dev
multiqc=bfq-dev
microbiome=bfq-dev
OPTIND=3
while getopts d:r:s:t:q:m: flag
do
    case "${flag}" in
        d) gcfdb=${OPTARG};;
        r) rna_seq=${OPTARG};;
        s) single_cell=${OPTARG};;
        t) gcf_tools=${OPTARG};;
        q) multiqc=${OPTARG};;
        m) microbiome=${OPTARG};;
    esac
done
if [ "$1" == "test" ];
then
	BUILD_ARGS="--build-arg GCFDB_BRANCH=$gcfdb --build-arg RNA_SEQ_BRANCH=$rna_seq --build-arg SINGLE_CELL_BRANCH=$single_cell --build-arg GCF_TOOLS_BRANCH=$gcf_tools --build-arg MULTIQC_BRANCH=$multiqc --build-arg MICROBIOME_BRANCH=$microbiome"
else
	BUILD_ARGS=""
fi

DOCKERFILE="dockerfile-$1"
echo "BUILDING WITH DOCKERFILE: $DOCKERFILE"

sudo docker build -t gcfntnu/bfq:$2 . -f $DOCKERFILE $BUILD_ARGS
sudo docker push gcfntnu/bfq:$2
