#!/bin/bash

INSTRUMENTS_DIR=/home/geiramh/instruments
OUTPUT_DIR=/home/geiramh/bfq/
INDEX_DIR=/home/geiramh/indexes/
MASKED_INDEX_DIR=/home/geiramh/Documents/masked-genomes
FLOWCELLMANAGER_DIR=/home/geiramh/flowcellmanager
BFQ_VERSION_DIR=/home/geiramh/bfq_version

mkdir -p $BFQ_VERSION
echo "gcfntnu/bfq:dev-test" > $BFQ_VERSION_DIR/bfq.version
sudo docker kill bfq
sudo docker build -t bfq-pipe .
sudo docker run -d --privileged --rm --name=bfq  -v $MASKED_INDEX_DIR:/masked-genomes/ -v $OUTPUT_DIR:/bfq -v $INSTRUMENTS_DIR:/instruments -v $INDEX_DIR:/data -v $FLOWCELLMANAGER_DIR:/flowcellmanager -v $BFQ_VERSION_DIR:/bfq_version bfq-pipe

