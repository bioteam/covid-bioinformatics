#!/bin/bash
# Run hmmscan against an HMM database with one or more sequences

seqs=FALSE
db=FALSE
UNDERSCORE='_'

while getopts ":s:d:" opt
   do
     case $opt in
        s ) seqs=$OPTARG;;
        d ) db=$OPTARG;;
     esac
done

if [[ $seqs == FALSE || $db == FALSE ]]
then
    echo "Usage: $0 -s '<sequence files>' -d <HMM database>"
    exit
fi

for seq in $seqs
do
    if [ ! -f $seq ]
    then
        echo Sequence file "$seq" not found
        exit
    fi
    seqname=$( basename $seq | cut -d '.' -f 1 )
    dbname=$( basename $db | cut -d '.' -f 1 )
    out=${seqname}${UNDERSCORE}${dbname}.hmmscan
    if [ ! -f $out ]
    then
        cmd="hmmscan -o $out $db $seq"
	    echo Command is $cmd
        $cmd
    fi
done
