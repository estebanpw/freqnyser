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

float difpor = 50;
uint64_t min_freq = 100;

std::list<Significative> * compare_kmers(uint64_t * x, uint64_t * y, uint64_t total);
bool compare_significative(const Significative& first, const Significative& second)
{
    return ((first.difpor) > (second.difpor));
}

int main(int argc, char ** av){

    uint64_t i, custom_kmer = 2;

    // @@@@@@@@@@@@@@@@@@@@@@@ Parameter init stuff

    if(argc != 6) { fprintf(stderr, "ERROR: Use %s noh-mers h-mers ksize out mindifpor\n", av[0]); exit(-1); }

    FILE * file_x = fopen(av[1], "rb");
    FILE * file_y = fopen(av[2], "rb");

    if(file_x == NULL) terror("No file X"); 
    if(file_y == NULL) terror("No file Y");
    custom_kmer = (uint64_t) atoi(av[3]);

    FILE * out = fopen(av[4], "wt");
    if(out == NULL) terror("No out file");

    difpor = atof(av[5]);

    uint64_t total = (uint64_t) 4 << (2*custom_kmer);

    uint64_t * mers_table_x = (uint64_t *) calloc(total, sizeof(uint64_t));
    if(mers_table_x == NULL) terror("Could not allocate first mers table");

    uint64_t * mers_table_y = (uint64_t *) calloc(total, sizeof(uint64_t));
    if(mers_table_y == NULL) terror("Could not allocate second mers table");
    
    for(i=0; i<custom_kmer; i++){
        total = (uint64_t) 4 << (2*i);
        fread(&mers_table_x[0], sizeof(uint64_t), total, file_x);
        fread(&mers_table_y[0], sizeof(uint64_t), total, file_y);

        std::list<Significative> * interesting = compare_kmers(mers_table_x, mers_table_y, total);
        interesting->sort(compare_significative);

        std::cout << "K=" << i+1 << "------------------------------\n";

        std::list<Significative>::iterator it; int j = 0;
        for (it=interesting->begin(); it!=interesting->end(); ++it){

            char word[custom_kmer+1]; 
            perfect_hash_to_word(word, (*it).hash, (uint32_t) i+1);
            word[i+1] = '\0';

            if(j>20) std::cout << (*it).hash << " " << word << " " << (*it).noh_freq << " " << (*it).h_freq << " " << (*it).tx << " " << (*it).ty << "\t\t" << (*it).difpor << "\n";
            fprintf(out, "%" PRIu64" %s %Le %Le %" PRIu64" %" PRIu64" %f\n", (*it).hash, word, (*it).noh_freq, (*it).h_freq, (*it).tx, (*it).ty, (*it).difpor);
            
            ++j;
        }
        std::cout << "Significative found " << interesting->size() << "\n";
        fprintf(out, "# %" PRIu64"\n", interesting->size());

        interesting->clear();

    }
        

    fclose(file_x);
    fclose(file_y);
    fclose(out);
    free(mers_table_x);
    free(mers_table_y);
    
    return 0;
}

std::list<Significative> * compare_kmers(uint64_t * x, uint64_t * y, uint64_t total){
    uint64_t i, c_hash = 0, sum_x = 0, sum_y = 0;

    std::list<Significative> * sign_kmers = new std::list<Significative>();

    for(i=0; i<total; i++) { sum_x += x[i]; sum_y += y[i]; }

    
    for(i=0; i<total; i++)
    {

        Significative s;
        s.hash = c_hash;
        s.noh_freq = (long double) x[i] / (long double) sum_x;
        s.h_freq =  (long double) y[i] / (long double) sum_y;
        s.difpor = (float) ((s.h_freq - s.noh_freq)/(s.noh_freq)*100);
        s.tx = x[i];
        s.ty = y[i];

        if(s.difpor > difpor && y[i] > min_freq)
        {
            
            sign_kmers->push_back(s);

        }


        ++c_hash;
    }

    return sign_kmers;
}