#!/bin/bash

# AUTHOR	Marcos SR
#
# USAGE 	./batch_apply_massgen_to_folder.sh <input_folder> <output_folder>
#
# DESCRIPTION	Given an input folder, (1)creates a new folder called massgen_<kmer_size> inside the <output_folder>, (2)executes
#		Esteban's massive-gen (generates the kmers frequency table of a fasta file) for each .fasta inside <input_folder> into massgen_<kmer_size>

input_folder=$1/*.fasta
output_folder=$2
kmer_size=5

massgen_folder=${output_folder}/massgen_${kmer_size}
mkdir -p $massgen_folder

echo "> Generating massgen files..."
for fasta in $input_folder
do
	echo "[INFO] $fasta"
	name=$(basename "$fasta")
	name=${name%.*}
	echo "[INFO] $name"
	echo "/ssd/estebanpw/irene/freqnyser/bin/massive-gen -query $fasta -kmer $kmer_size -out ${massgen_folder}/${name}.massgen"
	/ssd/estebanpw/irene/freqnyser/bin/massive-gen -query $fasta -kmer $kmer_size -out ${massgen_folder}/${name}.massgen
done
