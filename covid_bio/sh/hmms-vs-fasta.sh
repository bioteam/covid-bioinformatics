#!/bin/bash
# Run HMMs against a Fasta file using hmmsearch,
# create 'table' output files like "S-aa_nr.tblout".

fasta=FALSE
hmms=FALSE
evalue=FAwLSE
UNDERSCORE='_'

while getopts ":h:f:e:" opt
   do
     case $opt in
        h ) hmms=$OPTARG;;
        f ) fasta=$OPTARG;;
        e ) evalue=$OPTARG;;
     esac
done

if [[ $fasta == FALSE || $hmms == FALSE ]]
then
    echo "Usage: $0 -f <fasta file> -h '<HMM files>' [-e <e value>]"
    exit
fi

if [ ! -f $fasta ]
then
    echo Fasta file "$fasta" not found
    exit
fi

for hmm in $hmms
do
    if [ ! -f $hmm ]
    then
        echo HMM "$hmm" not found
        exit
    fi
    base=$( echo $hmm | cut -d '.' -f 1 )
    out=${base}${UNDERSCORE}${fasta}.tblout
    if [ ! -f $out ]
    then
        echo HMM: $hmm Fasta: $fasta
        if [[ $evalue == FALSE ]]
        then
            cmd="hmmsearch --noali --tblout $out $hmm $fasta"
        else
            cmd="hmmsearch --noali -E $evalue --tblout $out $hmm $fasta"
        fi
	    echo Command is $cmd
        $cmd
        echo Output is $out
    else
	echo File $out exists
    fi
done
