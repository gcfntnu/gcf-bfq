#!/bin/bash

if [ "$#" -ne 2 ] ; 
then
    echo -e "USAGE IS: \n\n\t\t ./build-tag-push.sh [prod|test|base] [docker-tag]\n\n"
	exit 1
fi

DOCKERFILE="dockerfile-$1"
echo "BUILDING WITH DOCKERFILE: $DOCKERFILE"

sudo docker build -t gcfntnu/bfq:"$2" . -f $DOCKERFILE
sudo docker push gcfntnu/bfq:"$2"
