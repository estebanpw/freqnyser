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
#define MAX_MOTIF 5000
#define BUFFER 1000
#define RANGE 2
#define MINIMUM_SIZE 50

float THRESHOLD = 0.8;
uint64_t min_chain_dist = 5;

uint64_t get_seq_len(FILE * f);
uint64_t load_seq(FILE * f, char * seq);
void reverse_complement(char * origin, char * dest);
void init_args(int argc, char ** av, FILE ** query, FILE ** out_results, char * motif, uint64_t * max_chain_dist);
void find_chains(FILE * query, uint64_t max_chain_dist, char * motif, FILE * out_re);
int search_chain(uint64_t * current_position, uint64_t * position, char * seq_x, uint64_t query_l, char * motif);

int main(int argc, char ** av){

    // @@@@@@@@@@@@@@@@@@@@@@@ Parameter init stuff

    FILE * query = NULL, * out_re = NULL;
    uint64_t max_chain_dist = 1000;
    char motif[MAX_MOTIF]; motif[0] = '\0';
    init_args(argc, av, &query, &out_re, motif, &max_chain_dist);
    min_chain_dist = 2*strlen(motif);

    fprintf(stdout, "Searching motifs [%s] of len [%d]!\n", motif, (int)strlen(motif));


	fprintf(out_re, "#KMER %s\n", motif);
    find_chains(query, max_chain_dist, motif, out_re);

    fclose(out_re);
    
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

void init_args(int argc, char ** av, FILE ** query, FILE ** out_results, char * motif, uint64_t * max_chain_dist){

    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           consecutive-finder -query [query] -out [out] -motif [motif]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -maxdist         [Integer:   k>1 (default 1500)]\n");
            fprintf(stdout, "           -pident          [Float:     0<f<1 (default 0.8)]\n");
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
            *out_results = fopen(av[pNum+1], "a");
            if(out_results==NULL) terror("Could not open output file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-motif") == 0){
            strcpy(motif, av[pNum+1]);
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-maxdist") == 0){
            *max_chain_dist = (uint64_t) atoi(av[pNum+1]);
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-pident") == 0){
            THRESHOLD = atof(av[pNum+1]);
            pNum+=2;
        }
        

        if(pNum < 1) ++pNum;
        
    }
    if(*query == NULL || *out_results == NULL || strlen(motif) == 0) terror("Need query, output and motif");
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

void find_chains(FILE * query, uint64_t max_chain_dist, char * motif, FILE * out_re){

    uint64_t query_l = get_seq_len(query);
    uint64_t current_position = 0, s1 = 0, s2 = 0, s3 = 0;
    uint64_t motif_len = strlen(motif);

    char * seq_x = NULL;
    if ((seq_x = (char *) calloc(query_l, sizeof(char))) == NULL) {
        terror("Could not allocate memory for sequence x");
    }

    load_seq(query, seq_x);


    while(current_position < query_l){

        //printf("Lets go %" PRIu64" @ %.*s\n", current_position, 10, &seq_x[current_position]);
        if(1 == search_chain(&current_position, &s1, seq_x, query_l, motif)){
            current_position = s1 + motif_len;
            //printf("Lets go 2\n");
            if(1 == search_chain(&current_position, &s2, seq_x, query_l, motif) && (s2-s1) < max_chain_dist && (s2-s1+motif_len) > min_chain_dist){
                current_position = s2 + motif_len;
                //printf("Lets go 3\n");
                if(1 == search_chain(&current_position, &s3, seq_x, query_l, motif) && (s3-s2) < max_chain_dist && (s3-s2+motif_len) > min_chain_dist){
                    

                    current_position = s2 + motif_len; // Go back to second chain to try from there
                    s3 = s3 + 2*motif_len; // Increase last so that it also takes extra DNA in the motif
                    if(s3 > query_l) s3 = query_l;

                    fprintf(stdout, "%" PRIu64 ", %" PRIu64 ", %" PRIu64 "\n", s1, s2, s3);
                    fprintf(out_re, "> %" PRIu64 ", %" PRIu64 ", %" PRIu64 "\n", s1, s2, s3);
                    fprintf(out_re, "%.*s\n", (int)(s3-s1+motif_len), &seq_x[s1] );

                }
            }
        }
        ++current_position;
    }

    free(seq_x);

}

int search_chain(uint64_t * current_position, uint64_t * position, char * seq_x, uint64_t query_l, char * motif){
    
    char c = '\n';
    uint64_t matches = 0, l = 0;
    uint64_t motif_pos = 0, motif_len = strlen(motif);
    float pident;



    while( *current_position < query_l ) {
        
        c = seq_x[(*current_position)++];
            
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') {

            // AAAAATTAGCCAA
            
            if(motif[motif_pos++] == c) ++matches;

            ++l;

            pident = (float) matches / (float) l;

            
            if(pident < THRESHOLD){
                matches = 0;
                motif_pos = 0;
                l = 0;
            }

            if(motif_pos == motif_len && pident >= THRESHOLD){
                
                *position = *current_position - motif_pos;
                return 1;
            }

        }else{ //It can be anything (including N, Y, X ...)            
            
            ++(*current_position);
            matches = 0;
            motif_pos = 0;
            l = 0;

        }
    }


    return 0;


}
