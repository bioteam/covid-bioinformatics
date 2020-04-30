#!/bin/bash
# Run HMMs against a Fasta file using hmmsearch,
# create 'table' output files like S-aa_E-aa.out. Usage:
# hmms-vs-fasta.sh -h 'NS*aa.hmm' -f nr

fasta=FALSE
hmms=FALSE
evalue=1.0e-10
UNDERSCORE='_'

while getopts ":h:f:" opt
   do
     case $opt in
        h ) hmms=$OPTARG;;
        f ) fasta=$OPTARG;;
     esac
done

if [[ $fasta == FALSE || $hmms == FALSE ]]
then
    echo "Usage: $0 -f <fasta file> -h '<HMM files>'"
    exit
fi

if [ ! -f $fasta ]; then
    echo Fasta file "$fasta" not found
    exit
fi

for hmm in $hmms
do
    if [ ! -f $hmm ]; then
        echo HMM "$hmm" not found
        exit
    fi
    echo HMM: $hmm Fasta:$fasta
    base=$( echo $hmm | cut -d '.' -f 1 )
    out="$base$UNDERSCORE$fasta.out"
    cmd="hmmsearch --noali -E $evalue --tblout $out $hmm $fasta"
    $cmd
    echo Output is $out
done
