#!/bin/bash

# AUTHOR	Marcos SR
#
# USAGE 	./batch_read_massgen_files.sh <input_folder> <output_folder>
#
# DESCRIPTION	Given an input folder, (1)replaces in all fasta files spaces by underscores, (2)using Estebans
#		massive-gen generates their kmers frequencies in a new folder, (3)then using Marcoss read-massgen-file
#		generates CSVs for each massgen file and each kmer size in the output folder in the newly created
#		chromosomes

input_folder=$1/*.fasta
output_folder=$2
kmer_size=5

echo "> Extracting each sequence from each multifasta"

simple_folder=${output_folder}/simple_fasta
mkdir -p $simple_folder

for multifasta in $input_folder
do
	echo "python /ssd/estebanpw/irene/freqnyser/src/multi2simple.py -i $multifasta -o $simple_folder"
	python /ssd/estebanpw/irene/freqnyser/src/multi2simple.py -i $multifasta -o $simple_folder
done

#echo "> Replacing all file spaces with underscores..."
# rename all with underscores instead of spaces
#for file in $input_folder
#do
#	mv "$file" ${file// /_}
#done

#primera parte --> separar
massgen_folder=${output_folder}/massgen_${kmer_size}
mkdir -p $massgen_folder

echo "> Generating massgen files..."
for fasta in $simple_folder/*.fasta
do
	echo "[INFO] $fasta"
	name=$(basename "$fasta")
	name=${name%.*}
	#name=${name// /_}
	echo "[INFO] $name"
	echo "/ssd/estebanpw/irene/freqnyser/bin/massive-gen -query $fasta -kmer $kmer_size -out ${massgen_folder}/${name}.massgen"
	/ssd/estebanpw/irene/freqnyser/bin/massive-gen -query $fasta -kmer $kmer_size -out ${massgen_folder}/${name}.massgen
done

# en output crear la carpeta done
done_folder=${output_folder}/chromosomes
mkdir -p $done_folder

echo "> Reading massgen files and generating CSVs..."
for massgen in ${massgen_folder}/*.massgen
do
	echo "[INFO] $massgen"
        name=$(basename $massgen)
        name=${name%.*}
        echo "[INFO] $name"

	for k in $(seq 2 $kmer_size)
	do
		echo "/ssd/estebanpw/irene/freqnyser/bin/read-massgen-file $massgen $k ${done_folder}/${name}_${k}.csv"
		chrom_folder=${done_folder}/chrom_${k}
		mkdir -p $chrom_folder
		/ssd/estebanpw/irene/freqnyser/bin/read-massgen-file $massgen $k ${chrom_folder}/${name}_${k}.csv
	done
done

tables_folder=${output_folder}/tables
mkdir -p $tables_folder

echo "python /ssd/estebanpw/irene/freqnyser/src/indexes2table.py -i ${done_folder}/chrom_2 -o $tables_folder"
python /ssd/estebanpw/irene/freqnyser/src/indexes2table.py -i ${done_folder}/chrom_2 -o $tables_folder

