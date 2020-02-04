if [ "$#" -eq 2 ] ; 
then
    MYTAG=$2
else
    MYTAG=0
fi

sudo docker build -t bfq-pipe . -f dockerfile-test --build-arg TAG=$MYTAG
#sudo docker tag bfq-pipe gcfntnu/bfq:"$@"
#sudo docker push gcfntnu/bfq:"$@"
sudo docker tag bfq-pipe gcfntnu/bfq:"$1"
sudo docker push gcfntnu/bfq:"$1"
