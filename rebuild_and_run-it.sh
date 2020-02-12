#!/bin/bash
INSTRUMENTS_DIR=/home/geiramh/instruments
INSTRUMENTS_DIR=instruments
OUTPUT_DIR=/home/geiramh/bfq/
INDEX_DIR=/home/geiramh/indexes/
CONFIG_DIR=/home/geiramh/config-bfq/
MASKED_INDEX_DIR=/home/geiramh/Documents/masked-genomes
FLOWCELLMANAGER_DIR=/home/geiramh/flowcellmanager
BFQ_VERSION_DIR=/home/geiramh/bfq_version
INSTR_ARCHIVE=/home/geiramh/archive-instruments

if [ "$#" -ne 1 ] ;
then
    echo -e "USAGE IS: \n\n\t\t ./rebuild_and_run-it.sh [prod|test]\n\n"
    exit 1
fi


DOCKERFILE="dockerfile-$1"
echo "BUILDING WITH DOCKERFILE: $DOCKERFILE"

mkdir -p $BFQ_VERSION_DIR
echo "gcfntnu/bfq:dev-test" > $BFQ_VERSION_DIR/bfq.version
sudo docker kill bfq
sudo docker build -t bfq-pipe . -f $DOCKERFILE
sudo docker run --privileged --rm -it --name=bfq-it  -v $CONFIG_DIR:/config -v $MASKED_INDEX_DIR:/masked-genomes/ -v $OUTPUT_DIR:/bfq -v $INSTRUMENTS_DIR:/instruments -v $INDEX_DIR:/data -v $FLOWCELLMANAGER_DIR:/flowcellmanager -v $BFQ_VERSION_DIR:/bfq_version -v $INSTR_ARCHIVE:/instrument-archive --expose 5000 --net="host" bfq-pipe /bin/bash
