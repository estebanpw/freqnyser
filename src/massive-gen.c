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
#include "structs.h"
#include "commonFunctions.h"
#define STARTING_SEQS 1000
#define VECTOR_L 2000
#define MAX_NAME 1000
#define BUFFER 1000
#define RANGE 2
#define MINIMUM_SIZE 50

// To reduce overhead in recursion
uint64_t query_l, db_l;
char * seq_x, * seq_y;
long double theta;
float THRESHOLD = 4;

uint64_t get_seq_len(FILE * f);
uint64_t load_seq(FILE * f, char * seq);
void reverse_complement(char * origin, char * dest);
void init_args(int argc, char ** av, FILE ** query, FILE ** out_results, uint64_t * custom_kmer);
void generate_distribution(uint64_t ** mers_table, uint64_t k, FILE * query, uint64_t query_l);


int main(int argc, char ** av){

    uint64_t custom_kmer = 2, query_l, i;

    // @@@@@@@@@@@@@@@@@@@@@@@ Parameter init stuff

    FILE * query = NULL, * out_results = NULL;
    init_args(argc, av, &query, &out_results, &custom_kmer);



    uint64_t ** mers_table = (uint64_t **) calloc(custom_kmer, sizeof(uint64_t **));
    if(mers_table == NULL) terror("Could not allocate mers table");

    for(i=0; i<custom_kmer; i++){
        uint64_t total = (uint64_t) 4 << (2*i);
        mers_table[i] = (uint64_t *) calloc(total, sizeof(uint64_t *));
        if(mers_table[i] == NULL) terror("Could not alocate second loop");
    }

    

    query_l = get_seq_len(query);

    fprintf(stdout, "[INFO] Generating frequency distributions\n");

    generate_distribution(mers_table, custom_kmer, query, query_l);

    fprintf(stdout, "\n[INFO] Completed\n");


    

    fclose(query);

    for(i=0; i<custom_kmer; i++){
        fwrite(mers_table[i], sizeof(uint64_t), (uint64_t)(4 << (2*i)), out_results);
        /*
        for(uint64_t j=0; j<(uint64_t)(4 << (2*i)); j++){
            printf("%" PRIu64"\n", mers_table[i][j]);
        }
        printf("----------for k=%" PRIu64"\n", i);
        */
        free(mers_table[i]);
    } 
    free(mers_table);
    fclose(out_results);
    
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

void init_args(int argc, char ** av, FILE ** query, FILE ** out_results, uint64_t * custom_kmer){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           massive-gen -query [query] -out [out]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -kmer       [Integer:   k>1 (default 2)]\n");
            fprintf(stdout, "           -lambda     [Float:     Scalar for minimum average threshold\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        else if(strcmp(av[pNum], "-query") == 0){
            *query = fopen(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open query file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-out") == 0){
            *out_results = fopen(av[pNum+1], "wt");
            if(out_results==NULL) terror("Could not open output file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-kmer") == 0){
            *custom_kmer = (uint64_t) atoi(av[pNum+1]);
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-lambda") == 0){
            THRESHOLD = atof(av[pNum+1]);
            pNum+=2;
        }
        

        if(pNum < 1) ++pNum;
        
    }
    if(*query == NULL || *out_results == NULL) terror("A query and output are required");
}

uint64_t get_seq_len(FILE * f) {
    char c = '\0';
    uint64_t l = 0;

    while(!feof(f)){
        c = getc(f);

        if(c == '>'){
            ++l;
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
            seq[l] = '>';
            ++l;
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

void generate_distribution(uint64_t ** mers_table, uint64_t k, FILE * query, uint64_t query_l){

    
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

    uint64_t masks[k];
    for(uint64_t i=0; i<k; i++){
        masks[i] = (1 << (2*(i+1))) - 1;
    }



    while( current_len < query_l ) {
        
        c = seq_x[current_len++];
            
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            
            if(word_size < k) ++word_size;

            if(c == 'A') { hash = (hash << 2) + 0; }
            if(c == 'C') { hash = (hash << 2) + 1; }
            if(c == 'G') { hash = (hash << 2) + 2; }
            if(c == 'T') { hash = (hash << 2) + 3; }

            
            for(uint64_t i=0; i<k; i++){
                
                if(word_size > i) ++mers_table[i][hash & masks[i]];

            }
            
            if(current_len % a_hundreth == 0) { 
                fprintf(stdout, "\r%" PRIu64"%%...", 1+100*current_len/query_l); 
                fflush(stdout);
            }

        }else{ //It can be anything (including N, Y, X ...)            
            
            hash = 0;
            word_size = 0;
            ++current_len;

        }
    }

    free(seq_x);

}