#!/bin/bash

if [ "$#" -lt 2 ] ; 
then
    echo  -e "Usage: \n" \
    "./build-tag-push.sh [prod|test|base] [docker-tag] [optional args]\n\n" \
    "Optional branch args for dockerfile-test:\n" \
    "-t gcf-tools \n" \
    "-q multiqc \n" \
    "-w gcf-workflows \n" 
	exit 1
fi

gcf_tools=bfq-dev
multiqc=bfq-dev
gcf_workflows=bfq-dev
OPTIND=3
while getopts t:q:w: flag
do
    case "${flag}" in
        t) gcf_tools=${OPTARG};;
        q) multiqc=${OPTARG};;
        w) gcf_workflows=${OPTARG};;
    esac
done
if [ "$1" == "test" ];
then
	BUILD_ARGS="--build-arg GCF_TOOLS_BRANCH=$gcf_tools --build-arg MULTIQC_BRANCH=$multiqc --build-arg GCF_WORKFLOWS_BRANCH=$gcf_workflows"
else
	BUILD_ARGS=""
fi

DOCKERFILE="dockerfile-$1"
echo "BUILDING WITH DOCKERFILE: $DOCKERFILE"

echo "sudo docker build -t gcfntnu/bfq:$2 . -f $DOCKERFILE $BUILD_ARGS"
sudo docker build -t gcfntnu/bfq:$2 . -f $DOCKERFILE $BUILD_ARGS
sudo docker push gcfntnu/bfq:$2
