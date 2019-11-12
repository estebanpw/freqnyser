
MOTIF=$1
mkdir MEGA

for i in /ssd/estebanpw/irene/freqnyser/src/fastahotspots/*.fasta 
do
	awk 'BEGIN{x=0; header=0; current=0;}{if(x==0) header=$0; if(x==1){ name="MEGA/file"current; printf("%s\n%s\n", header, $0) > name;} x++; current++; if(x==2){ x=0;} }' $i
	for j in MEGA/*
	do
		/ssd/estebanpw/irene/freqnyser/bin/one-motif-finder -query $j -out $j.motif -motif $MOTIF -pident 0.75
	done
	grep ">" MEGA/*.motif | awk '{print $1}' | uniq > $(basename $i)_splitted.motifs
	rm MEGA/*
done

rm -rf MEGA
