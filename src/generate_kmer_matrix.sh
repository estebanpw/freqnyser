#!/bin/bash

# AUTHOR	Marcos SR
#
# USAGE 	./generate_kmer_matrix.sh  <input_folder> <output_folder>
#
# DESCRIPTION	(1)Reads .massgen files from <input_folder> with read-massgen-file and writes the kmer frequency tables of <kmer_size>
#		in a .csv inside <output_folder> (2)to create the frequency table of each previous file with indexes2table.py

input_folder=$1
output_folder=$2
kmer_size=5

# en output crear la carpeta done
done_folder=${output_folder}/chromosomes
mkdir -p $done_folder

echo "> Reading massgen files and generating CSVs..."
for massgen in ${input_folder}/*.massgen
do
	echo "[INFO] $massgen"
        name=$(basename $massgen)
        name=${name%.*}
        echo "[INFO] $name"

	#for k in $(seq 2 $kmer_size)
	#do
		echo "/ssd/estebanpw/irene/freqnyser/bin/read-massgen-file $massgen $kmer_size ${done_folder}/${name}_${kmer_size}.csv"
		chrom_folder=${done_folder}/chrom_${kmer_size}
		mkdir -p $chrom_folder
		/ssd/estebanpw/irene/freqnyser/bin/read-massgen-file $massgen $kmer_size ${chrom_folder}/${name}_${kmer_size}.csv
	#done
done

#tables_folder=${output_folder}/tables
#mkdir -p $tables_folder
echo "> Generating table from frequency tables..."
echo "python3 /ssd/estebanpw/irene/freqnyser/src/indexes2table.py -i ${done_folder}/chrom_${kmer_size} -o ${output_folder}"
python3 /ssd/estebanpw/irene/freqnyser/src/indexes2table.py -i ${done_folder}/chrom_${kmer_size} -o ${output_folder}
