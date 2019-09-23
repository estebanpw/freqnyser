/*********

File        unikmers.c
Author      EPW <estebanpw@uma.es>
Description Computes HSPs using unique hits with a variable degree of uniqueness

USAGE       Usage is described by calling ./unikmers --help



**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
//#include <iostream>
//#include <list>
#include "structs.h"
#include "commonFunctions.h"
#define STARTING_SEQS 1000
#define VECTOR_L 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes
#define RANGE 2

// To reduce overhead in recursion
uint64_t query_l, db_l;
char * seq_x, * seq_y;


uint64_t get_seq_len(FILE * f);
uint64_t load_seq(FILE * f, char * seq);
void reverse_complement(char * origin, char * dest);
void init_args(int argc, char ** av, FILE ** query, FILE ** out_results, uint64_t * custom_kmer, uint64_t * selected_hash);
void generate_distribution(uint64_t * frequency_table, Dictionary * d, uint64_t k, FILE * query, uint64_t query_l);
void generate_plotting_vector(uint64_t * frequency_table, Dictionary * d, uint64_t k, uint64_t query_l, uint64_t selected_hash, FILE * out_results, uint64_t factor_scale);
void print_distribution(uint64_t * frequency_table, uint64_t k);


int main(int argc, char ** av){

    uint64_t custom_kmer = 2, query_l, selected_hash = 0;

    

    FILE * query = NULL, * out_results = NULL;
    init_args(argc, av, &query, &out_results, &custom_kmer, &selected_hash);

    uint64_t * frequency_table = (uint64_t *) calloc(1 << (2*custom_kmer), sizeof(uint64_t));
    if(frequency_table == NULL) terror("Could not allocate frequency table");

    query_l = get_seq_len(query);

    Dictionary * d = (Dictionary *) calloc(query_l, sizeof(Dictionary));
    if(d == NULL) terror("Could not generate dictionary, possibly because there is not enough RAM (16 bytes required per nucleotide)");

    fprintf(stdout, "[INFO] Using k size %" PRIu64", hash selected %" PRIu64"\n", custom_kmer, selected_hash);

    fprintf(stdout, "[INFO] Generating frequency distributions.\n");

    generate_distribution(frequency_table, d, custom_kmer, query, query_l);

    fprintf(stdout, "\n[INFO] Completed.\n");

    //if(custom_kmer <= 8) print_distribution(frequency_table, custom_kmer);

    generate_plotting_vector(frequency_table, d, custom_kmer, query_l, selected_hash, out_results, query_l / VECTOR_L);

    fclose(query);
    if(out_results != NULL && out_results != stdout) fclose(out_results);

    free(frequency_table);
    free(d);


    

    
    /*
    
    // Get file length
    
    db_l = get_seq_len(database);
    a_hundreth = MAX(1, (db_l/100));

    seq_y = NULL;
    if ((seq_y = (char *) calloc(db_l, sizeof(char))) == NULL) {
        terror("Could not allocate memory for sequence y");
    }

    load_seq(database, seq_y);

    fprintf(stdout, "[INFO] Computing database.\n");


    current_len = 0;

    curr_kmer[0] = reverse_kmer[0] = '\0';
    word_size = 0, word_size_rev = custom_kmer-1;

    
    while( current_len < db_l ) {
        
        c = seq_y[current_len++];

        
            
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            
            curr_kmer[word_size] = c;
            switch(c){
                case ('A'): reverse_kmer[word_size_rev] = (unsigned)'T';
                break;
                case ('C'): reverse_kmer[word_size_rev] = (unsigned)'G';
                break;
                case ('G'): reverse_kmer[word_size_rev] = (unsigned)'C';
                break;
                case ('T'): reverse_kmer[word_size_rev] = (unsigned)'A';
                break;
            }
            if(word_size_rev != 0) --word_size_rev;

            if(word_size < custom_kmer) ++word_size;
            
            if(current_len % a_hundreth == 0) { 
                fprintf(stdout, "\r%" PRIu64"%%...", 1+100*current_len/db_l); 
                fflush(stdout);
            }

        }else{ //It can be anything (including N, Y, X ...)

            if(c != '\n' && c != '>') {
                
                word_size = 0;
                word_size_rev = custom_kmer-1;
                ++current_len;

            } 
        }

        if(word_size == custom_kmer){
            
            //fprintf(stdout, "\t Y-F %.*s -> %"PRIu64" @p:%"PRIu64"\n", 32, curr_kmer, hashOfWord(&curr_kmer[0], custom_kmer, 0), current_len);
            //fprintf(stdout, "\t Y-R %.*s -> %"PRIu64" @p:%"PRIu64"\n", 32, reverse_kmer, hashOfWord(&reverse_kmer[0], custom_kmer, 0), current_len);
            //getchar();

            root = insert_AVLTree_y(root, hashOfWord(&curr_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used, FORWARD);
            root = insert_AVLTree_y(root, hashOfWord(&reverse_kmer[0], custom_kmer, 0), mp_AVL, &n_pools_used_AVL, current_len, mp, &n_pools_used, REVERSE);

            // Overlapping
            
            memmove(&curr_kmer[0], &curr_kmer[1], custom_kmer-1);
            memmove(&reverse_kmer[1], &reverse_kmer[0], custom_kmer-1);
            --word_size;
        }
    }

    fprintf(stdout, "[INFO] Database computed of length %" PRIu64".\n", current_len);
    fclose(database);


    */


    return 0;
}


void reverse_complement(char * origin, char * dest){
    int i;
    for(i=0; i<32; i++){
        if(origin[i] == 'A') dest[32-(i+1)] = 'T';
        if(origin[i] == 'C') dest[32-(i+1)] = 'G';
        if(origin[i] == 'G') dest[32-(i+1)] = 'C';
        if(origin[i] == 'T') dest[32-(i+1)] = 'A';
    }
}

void init_args(int argc, char ** av, FILE ** query, FILE ** out_results, uint64_t * custom_kmer, uint64_t * selected_hash){

    int pNum = 0;
    uint64_t word_l = 0;
    *out_results = stdout;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           freqgen -query [query]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -out        [outfile] (default stdout)\n");
            fprintf(stdout, "           -kmer       [Integer:   k>1 (default 2)]\n");
            fprintf(stdout, "           -word       [String: word to plot (default AA)\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        else if(strcmp(av[pNum], "-query") == 0){
            *query = fopen64(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open query file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-out") == 0){
            *out_results = fopen(av[pNum+1], "wt");
            if(out_results==NULL) terror("Could not open output database file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-kmer") == 0){
            *custom_kmer = (uint64_t) atoi(av[pNum+1]);
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-word") == 0){

            uint64_t j, hash = 0;
            word_l = strlen(av[pNum+1]);
            for(j=0; j<strlen(av[pNum+1]); j++){

                if(av[pNum+1][j] == 'A') hash = (hash << 2) + 0;
                if(av[pNum+1][j] == 'C') hash = (hash << 2) + 1;
                if(av[pNum+1][j] == 'G') hash = (hash << 2) + 2;
                if(av[pNum+1][j] == 'T') hash = (hash << 2) + 3;
            }
            
            *selected_hash = hash;
            
            pNum+=2;
        }
        

        if(pNum < 1) ++pNum;
        
    }
    
    if(word_l != *custom_kmer) terror("Kmer size and word size must be equal. Use -kmer <n> and -word <x1x2..xN>");
    if(*query==NULL) terror("A query is required");
}

uint64_t get_seq_len(FILE * f) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            ++l;
        }
    }


    rewind(f);
    return l;
}

uint64_t load_seq(FILE * f, char * seq) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            while(c != '\n') c = getc(f);
        }
        c = toupper(c);
        if(c >= 'A' && c <= 'Z'){
            seq[l] = c;
            ++l;
        }
    }


    rewind(f);
    return l;
}

void generate_distribution(uint64_t * frequency_table, Dictionary * d, uint64_t k, FILE * query, uint64_t query_l){

    
    //char curr_kmer[k];
    //curr_kmer[0] = '\0';
    
    uint64_t word_size = 0;

    //To hold all information related to database
    uint64_t current_len = 0;   

    // Variables to read kmers
    char c = 'N'; //Char to read character
    
    // Read full sequence
    uint64_t a_hundreth = MAX(1, (query_l/100));

    seq_x = NULL;
    if ((seq_x = (char *) calloc(query_l, sizeof(char))) == NULL) {
        terror("Could not allocate memory for sequence x");
    }

    load_seq(query, seq_x);

    uint64_t hash = 0;

    uint64_t mask = (1 << (2*k)) - 1;



    while( current_len < query_l ) {
        
        c = seq_x[current_len++];
            
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            
            //curr_kmer[word_size] = c;
            if(word_size < k) ++word_size;

            if(c == 'A') hash = (hash << 2) + 0;
            if(c == 'C') hash = (hash << 2) + 1;
            if(c == 'G') hash = (hash << 2) + 2;
            if(c == 'T') hash = (hash << 2) + 3;

            hash = hash & mask;
            
            if(current_len % a_hundreth == 0) { 
                fprintf(stdout, "\r%" PRIu64"%%...", 1+100*current_len/query_l); 
                fflush(stdout);
            }

        }else{ //It can be anything (including N, Y, X ...)

            if(c != '\n' && c != '>') {
                
                hash = 0;
                word_size = 0;
                ++current_len;

            } 
        }

        if(word_size == k){
            
            //fprintf(stdout, "\t X-F %.*s -> %" PRIu64" @p:%" PRIu64"\n", (int) k, curr_kmer, hash, current_len);
            //getchar();

            //if(hash >= (uint64_t) (1 << (2*k))) fprintf(stdout, "REPORTED!!!!!!!\n");
            ++frequency_table[hash]; // True distribution
            
            d[current_len - k].hash = hash;
            d[current_len - k].position = current_len - k;

            // Overlapping
            //memmove(&curr_kmer[0], &curr_kmer[1], k-1);
            //--word_size;
        }
    }

    free(seq_x);

}

void generate_plotting_vector(uint64_t * frequency_table, Dictionary * d, uint64_t k, uint64_t query_l, uint64_t selected_hash, FILE * out_results, uint64_t factor_scale){

    // Format is:
    // frequency of the desired string
    // length of the sequence
    // kmers @ pixel 1
    // kmers @ pixel 2
    // etc.

    fprintf(out_results, "%" PRIu64"\n", frequency_table[selected_hash]);
    fprintf(out_results, "%" PRIu64"\n", query_l);

    uint64_t i, j = 0, current;
    for(i=0; i<VECTOR_L; i++){
        
        current = 0;

        while(d[j].position < (i+1)*factor_scale && j < query_l){

            if(d[j].hash == selected_hash){
                ++current;
            }

            ++j;
            
        }
        fprintf(out_results, "%" PRIu64"\n", current);

    }
}

void print_distribution(uint64_t * frequency_table, uint64_t k){
    uint64_t i;
    char word[k+1];
    
    for( i=0; i< (uint64_t) (1 << (2*k)); i++ ){

        perfect_hash_to_word(word, i, k);
        word[k] = '\0';
        fprintf(stdout, "%s : (%" PRIu64") \t\t %" PRIu64"\n", word, i, frequency_table[i]);
    }
}