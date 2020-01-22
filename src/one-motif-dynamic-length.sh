#!/bin/bash

# use should be ./this fastain folderoutput 5 20 A
# This will look for all poly A from size 5 to 20 


INPUT=$1
OUTPUT=$2
BEGIN=$3
END=$4
MOTIFBASE=$5


name=$(basename $INPUT)

for (( i=$BEGIN; i<=$END; i++  ))
do
	mymotif=""
	for (( j=0; j<$i; j++ ))
	do
		mymotif=${mymotif}$MOTIFBASE
	done

	echo "$mymotif"
	/ssd/estebanpw/irene/freqnyser/bin/one-motif-finder -query $INPUT -out $OUTPUT/poly${MOTIFBASE}-$i-$name -motif $mymotif -pident 0.9

done
