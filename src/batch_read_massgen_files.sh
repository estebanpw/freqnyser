
#./read-massgen-file /home/marcossr/hotspots/kmers_freq_in_hs/chr1.massgen 3 /home/estebanpw/asd.txt


# generar los massgen de cada chr1.fasta con massive-gen
# generar las frec abs y relativas con read_massgen_files

input_folder=/ssd/estebanpw/irene/freqnyser/src/fastahotspots/*.fasta
kmer_size=5
output_folder=/home/marcossr/hotspots/kmers_freq_in_hs

for fasta in $input_folder
do
	echo "$fasta"
	name=$(basename $fasta)
	name=${name%.*}
	#echo "$name"
	echo "/ssd/estebanpw/irene/freqnyser/bin/massive-gen -query $fasta -kmer $kmer_size -out ${output_folder}/${name}_${kmer_size}.massgen"
	/ssd/estebanpw/irene/freqnyser/bin/massive-gen -query $fasta -kmer $kmer_size -out ${output_folder}/${name}_${kmer_size}.massgen
done

# en output crear la carpeta done
done_folder=${output_folder}/done
mkdir $done_folder

for massgen in $output_folder/*.massgen
do
	echo "$massgen"
        name=$(basename $massgen)
        name=${name%.*}
        #echo "$name"

	for k in $(seq 2 $kmer_size)
	do
		echo "/ssd/estebanpw/irene/freqnyser/bin/read-massgen-file $massgen $k ${done_folder}/${name}_${k}.csv"
		/ssd/estebanpw/irene/freqnyser/bin/read-massgen-file $massgen $k ${done_folder}/${name}_${k}.csv
	done
done
