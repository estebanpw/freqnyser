#!/bin/bash

if [ $# != 3 ]; then
    echo "***ERROR*** Use: $0 chromo hotspots k"
    exit -1
fi


chromo=$1
hotspots=$2
k=$3

chromoName=$(basename "$chromo")
chromoName="${chromoName%.*}"
hotspotName=$(basename "$hotspots")
hotspotName="${hotspotName%.*}"

echo "$chromoName $hotspotName"

../bin/massive-gen -query $chromo -kmer $k -out $chromoName-$k.kmers 
../bin/massive-gen -query $hotspots -kmer $k -out $hotspotName-$k.kmers 


../bin/compare-gen $chromoName-$k.kmers $hotspotName-$k.kmers $k $chromoName-$hotspotName-$k.comparison 50


rm $chromoName-$k.kmers
rm $hotspotName-$k.kmers
