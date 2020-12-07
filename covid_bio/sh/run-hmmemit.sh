#!/bin/bash
# Run hmmemit with one or more input HMMs

hmms=FALSE
UNDERSCORE='_'

while getopts ":h:" opt
   do
     case $opt in
        h ) hmms=$OPTARG;;
     esac
done

if [[ $hmms == FALSE ]]
then
    echo "Usage: $0 -h '<HMM files>'"
    exit
fi

for hmm in $hmms
do
    if [ ! -f $hmm ]
    then
        echo HMM "$hmm" not found
        exit
    fi
    base=$( basename $hmm | cut -d '.' -f 1 )
    out=${base}${UNDERSCORE}hmmemit.fa
    if [ ! -f $out ]
    then
        cmd="hmmemit -c -o $out $hmm"
	    echo Command is "$cmd"
        $cmd
    fi
done
