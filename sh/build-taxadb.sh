#!/bin/bash
# taxadb - https://github.com/HadrienG/taxadb
# pip3 install --user taxadb
# export PATH=~/.local/bin:$PATH
# 

TAXHOME=${HOME}/taxonomy

cd $TAXHOME

if [[ -d taxadb ]]
then
    mv -f taxadb taxadb.bak
fi

for f in taxadb.sqlite taxadb.sqlite-shm taxadb.sqlite-wal
do
    if [[ -f $f ]]
    then
        mv -f $f $f.bak
    fi
done

taxadb download --outdir taxadb --type full --quiet
taxadb create --input taxadb --dbname taxadb.sqlite --chunk 100 --division full --fast --quiet

rm -fr *.bak
