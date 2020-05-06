#!/bin/bash
# taxadb - https://github.com/HadrienG/taxadb
# pip3 install --user taxadb
# export PATH=~/.local/bin:$PATH
#

TAXHOME=${HOME}/taxonomy

cd $TAXHOME

if [[ -d taxadb ]]
then
    mv taxadb taxadb.bak
    mkdir taxadb
fi

taxadb download --outdir taxadb --type prot
taxadb create --input taxadb --dbname taxadb.sqlite --chunk 100 --division prot

rm -fr taxadb.bak
