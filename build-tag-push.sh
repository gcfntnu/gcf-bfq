sudo docker build -t bfq-pipe .
sudo docker tag bfq-pipe gcfntnu/bfq:"$@"
sudo docker push gcfntnu/bfq:"$@"
