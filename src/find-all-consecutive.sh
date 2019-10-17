#!/bin/bash
if [ $# != 2 ]; then
    echo "***ERROR*** Use: $0 <seq.fasta> <ksize>"
    exit -1
fi

SEQUENCE=$1
KSIZE=$2
location=/ssd/estebanpw/irene/freqnyser/bin/consecutive-finder
name=$(basename $SEQUENCE)

grep -v ">" $SEQUENCE | tr -d '\n' | awk -v name="$name" -v seq="$SEQUENCE" -v path="$location" -v k="$KSIZE" '

{
	
	hyper=$1;
	mylen=length(hyper);
	for (i=0; i<mylen; i++){

		kmer = substr(hyper,i,k);
		if(checked[kmer] == ""){
			checked[kmer] = kmer;
			system(path " -query " seq " -out " name "-all-motifs.txt -motif " kmer)
			getline val < "-"
		}	
	}
}
'

