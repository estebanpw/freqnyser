
MOTIF=$1

mkdir motif$MOTIF

cd motif$MOTIF
 
mkdir splitted
mkdir notsplitted

cd notsplitted

for i in /ssd/estebanpw/irene/freqnyser/src/fastahotspots/*.fasta ; do /ssd/estebanpw/irene/freqnyser/bin/one-motif-finder -query $i -out notsplitted_$(basename $i).motifs -motif $MOTIF -pident 0.75 ; done

cd ..

for i in /home/estebanpw/data/chromosomes/homo_sapiens/HOMSA.Chr.* ; do /ssd/estebanpw/irene/freqnyser/bin/one-motif-finder -query $i -out fulldna_$(basename $i).motifs -motif $MOTIF -pident 0.75 ; done

cd splitted

cp /ssd/estebanpw/irene/freqnyser/src/find-motif-in-splitted-hotspots.sh .

echo "[RUN] ./find-motif-in-splitted-hotspots.sh $MOTIF"
./find-motif-in-splitted-hotspots.sh $MOTIF

cd ..
