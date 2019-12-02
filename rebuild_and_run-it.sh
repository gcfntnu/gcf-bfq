#INSTRUMENTS_DIR=/home/geiramh/instruments
INSTRUMENTS_DIR=instruments
OUTPUT_DIR=/home/geiramh/bfq/
INDEX_DIR=/home/geiramh/indexes/
CONFIG_DIR=/home/geiramh/config-bfq/
MASKED_INDEX_DIR=/home/geiramh/Documents/masked-genomes
FLOWCELLMANAGER_DIR=/home/geiramh/flowcellmanager
BFQ_VERSION_DIR=/home/geiramh/bfq_version
INSTR_ARCHIVE=/home/geiramh/archive-instruments

mkdir -p $BFQ_VERSION_DIR
echo "gcfntnu/bfq:dev-test" > $BFQ_VERSION_DIR/bfq.version
sudo docker kill bfq
sudo docker build -t bfq-pipe .
sudo docker run --rm -it --name=bfq-it  -v $CONFIG_DIR:/config -v $MASKED_INDEX_DIR:/masked-genomes/ -v $OUTPUT_DIR:/bfq -v $INSTRUMENTS_DIR:/instruments -v $INDEX_DIR:/data -v $FLOWCELLMANAGER_DIR:/flowcellmanager -v $BFQ_VERSION_DIR:/bfq_version -v $INSTR_ARCHIVE:/instrument-archive --expose 5000 --net="host" bfq-pipe /bin/bash
