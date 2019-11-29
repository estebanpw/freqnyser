/*********

AUTHOR		Marcos SR

USAGE		./read-massgen-file <input file> <custom kmer> <out file>
DESCRIPTION	Read a massive-gen generated file with parameter 12 and generates a custom kmer frequency table

**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <algorithm>
#include <set>
#include <list>
#include "structs.h"
#include "commonFunctions.h"
#define STARTING_SEQS 1000
#define VECTOR_L 2000
#define MAX_NAME 1000
#define BUFFER 1000
#define RANGE 2
#define MINIMUM_SIZE 50


int main(int argc, char ** av){

    uint64_t i, custom_kmer = 2, total_amount = 0, seq_length = 0;

    if(argc != 4) { fprintf(stderr, "ERROR: Use %s massgen_file ksize out_file\n", av[0]); exit(-1); }

    FILE * massgen_file = fopen(av[1], "rb");
    if(massgen_file == NULL) terror("No file");

    custom_kmer = (uint64_t) atoi(av[2]);

    FILE * out_file = fopen(av[3], "wt");
    if(out_file == NULL) terror("No out file");

    uint64_t total_ksize = (uint64_t) 4 << (2*custom_kmer);

    uint64_t * kmers_freq_abs = (uint64_t *) calloc(total_ksize, sizeof(uint64_t));
    if(kmers_freq_abs == NULL) terror("Could not allocate first kmers absolute table");

    for(i=0; i<custom_kmer; i++){
        total_ksize = (uint64_t) 4 << (2*i);
        fread(&kmers_freq_abs[0], sizeof(uint64_t), total_ksize, massgen_file);
        //std::cout << "K=" << i+1 << "------------------------------\n";
    }

    char word[custom_kmer+1];
    word[custom_kmer] = '\0';
    for(uint64_t i=0; i<total_ksize; i++){
    	total_amount += kmers_freq_abs[i];
    }

    
    // LEER LA MIERDA DE LA LONGITUD DE LA SECUENCIA QUE ESTA AL FINAL DEL FICHERO, JODER
    while(!feof(massgen_file)) fread(&seq_length, sizeof(uint64_t), 1, massgen_file);
    std::cout << "[LENGTH] " <<seq_length << "\n";
    
    fprintf(out_file, "kmer,abs_freq,rel_freq\n");

    // !!!!!!!!!!!PROBABLY NOT THE BEST WAY
    fprintf(out_file, "length,%" PRIu64",%" PRIu64"\n", seq_length, seq_length);


    //fprintf(out_file, "kmer,abs_freq,rel_freq\n");

    for(uint64_t i=0; i<total_ksize; i++){
    	perfect_hash_to_word(word, i, (uint32_t) custom_kmer);
    	//std::cout << word << ", " << kmers_freq_abs[i] << ", " << ((float)kmers_freq_abs[i]/(float)total_amount) << "\n";
    	fprintf(out_file, "%s,%" PRIu64",%f\n", word, kmers_freq_abs[i], ((float)kmers_freq_abs[i]/(float)total_amount));
	//fprintf(out_file, "%s%, f\n", word, ((float)kmers_freq_abs[i]/(float)total_amount));
    }

    fclose(massgen_file);
    fclose(out_file);
    free(kmers_freq_abs);

    return 0;
}
