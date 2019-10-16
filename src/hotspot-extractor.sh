#!/bin/bash

if [ $# != 2 ]; then
    echo "***ERROR*** Use: $0 <chromo> <hotspots file>"
    exit -1
fi


chr=$1
chrhotspot=$2
chrname=$(basename "$chrhotspot")

awk '{print $2,$3}' $chrhotspot > lines



while IFS= read -r line; do

	values=(${line// / })

	echo ">${values[0]} to ${values[1]}" >> allhotspots-$chrname.fasta
	./region-DNA-extract.sh ${values[0]} ${values[1]} $chr >> allhotspots-$chrname.fasta
	

done < lines




rm lines
