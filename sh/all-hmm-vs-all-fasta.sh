#!/bin/bash
#
# Run all HMMs against all Fasta files using hmmsearch,
# create 'table' output files like S-aa_E-aa.out. Usage:
# all-hmm-vs-all-fasta.sh -h 'NS*aa.hmm' -f 'NS*aa.fasta'

evalue=1.0e-10
UNDERSCORE='_'

while getopts ":h:f:" opt
   do
     case $opt in
        h ) hmms=$OPTARG;;
        f ) fastas=$OPTARG;;
     esac
done

for hmm in $hmms
do
    if [ ! -f $hmm ]; then
        echo HMM "$hmm" not found
        exit
    fi
    echo HMM is $hmm
    hfile=$(basename $hmm)
    hbase=$(echo $hfile | cut -d '.' -f 1 )
    for fasta in $fastas
    do
        if [ ! -f $fasta ]; then
            echo Fasta "$fasta" not found
            exit
        fi
        echo Fasta is $fasta
        ffile=$(basename $fasta)
        fbase=$(echo $ffile | cut -d '.' -f 1 )
        out="${hbase}${UNDERSCORE}${fbase}.out"
        cmd="hmmsearch --noali -E $evalue --tblout $out $hmm $fasta"
        $cmd
        echo Output is $out
    done
done