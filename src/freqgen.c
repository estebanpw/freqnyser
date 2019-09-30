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
#include <list>
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
std::list<Hotspot> * load_hotspots(FILE * hotspots);
void reverse_complement(char * origin, char * dest);
char * init_args(int argc, char ** av, FILE ** query, FILE ** hotspots, FILE ** out_results, uint64_t * custom_kmer, uint64_t * selected_hash);
void generate_distribution(uint64_t * frequency_table, uint64_t * frequency_nucl, Dictionary * d, uint64_t k, FILE * query, uint64_t query_l);
void generate_plotting_vector(uint64_t * frequency_table, uint64_t * frequency_nucl, Dictionary * d, uint64_t k, uint64_t query_l, uint64_t selected_hash, FILE * out_results, uint64_t factor_scale);
std::list<Region> * detect_density_regions(Dictionary * d, uint64_t * frequency_table, uint64_t query_l, uint64_t selected_hash, FILE * out_results);
void check_density_hotspots(Dictionary * d, uint64_t * frequency_table, std::list<Hotspot> * hotspot_list, uint64_t query_l, uint64_t selected_hash);
bool compare_regions(const Region& first, const Region& second);
bool compare_hotspots(const Hotspot& first, const Hotspot& second);
void write_regions(Dictionary * d, uint64_t query_l, uint64_t selected_hash, FILE * out_results);
void write_average_regions(Dictionary * d, uint64_t query_l, uint64_t selected_hash, FILE * out_results);
void print_distribution(uint64_t * frequency_table, uint64_t k);


int main(int argc, char ** av){

    uint64_t custom_kmer = 2, query_l, selected_hash = 0;    

    // @@@@@@@@@@@@@@@@@@@@@@@ Parameter init stuff

    FILE * query = NULL, * out_results = NULL, * raw_out = NULL, * density_out = NULL, * average_out = NULL, * hotspots = NULL;
    char * exp_name = init_args(argc, av, &query, &hotspots, &out_results, &custom_kmer, &selected_hash);
    char out_name[MAX_NAME]; out_name[0] = '\0'; 



    uint64_t * frequency_table = (uint64_t *) calloc(1 << (2*custom_kmer), sizeof(uint64_t));
    if(frequency_table == NULL) terror("Could not allocate frequency table");
    uint64_t frequency_nucl[4] = {0, 0, 0, 0};

    query_l = get_seq_len(query);

    Dictionary * d = (Dictionary *) calloc(query_l, sizeof(Dictionary));
    if(d == NULL) terror("Could not generate dictionary, possibly because there is not enough RAM (16 bytes required per nucleotide)");
    uint64_t i; for(i=0; i<query_l; i++){ d[i].position = -1; } // Initialize

    fprintf(stdout, "[INFO] Using k size %" PRIu64", hash selected %" PRIu64"\n", custom_kmer, selected_hash);

    fprintf(stdout, "[INFO] Generating frequency distributions.\n");

    generate_distribution(frequency_table, frequency_nucl, d, custom_kmer, query, query_l);

    fprintf(stdout, "\n[INFO] Completed.\n");

    
    // @@@@@@@@@@@@@@@@@@@@@@@ Automatic density detector
    if(hotspots == NULL){

        strcpy(&out_name[0], exp_name); strcat(out_name, "_raw.vec"); raw_out = fopen(out_name, "wt"); if(raw_out == NULL) terror("Could not open output vector file (1)");
        strcpy(&out_name[0], exp_name); strcat(out_name, "_density.vec"); density_out = fopen(out_name, "wt"); if(density_out == NULL) terror("Could not open output vector file (2)");
        strcpy(&out_name[0], exp_name); strcat(out_name, "_average.vec");average_out = fopen(out_name, "wt"); if(average_out == NULL) terror("Could not open output vector file (3)");

        generate_plotting_vector(frequency_table, frequency_nucl, d, custom_kmer, query_l, selected_hash, out_results, (uint64_t) ((float)query_l / (float) VECTOR_L) );

        std::list<Region> * density_list = detect_density_regions(d, frequency_table, query_l, selected_hash, out_results);

        fprintf(stdout, "[INFO] Sorting\n");

        density_list->sort(compare_regions);

        std::list<Region>::iterator it;
        
        fprintf(stdout, "[INFO] Top 50 results:\n");

        i = 0;
        fprintf(out_results, "EOP\n");
        for (it=density_list->begin(); it!=density_list->end(); ++it) {

            if((*it).p2 - (*it).p1 > MINIMUM_SIZE){
                float norm_score = (float) THRESHOLD * (float) theta + (*it).score; 
                float ratio = 0;
                if(theta != 0) ratio = norm_score / (float) theta;
                if(i<50 && hotspots != NULL) std::cout << " R: " << (*it).p1 << ", " << (*it).p2 << " - L " << (*it).p2 - (*it).p1 << " - S (norm) " << norm_score << " Ratio " << ratio << " Periods " << 1/norm_score << "/" << 1/theta << "\n";
                fprintf(density_out, "%" PRIu64",%" PRIu64",%" PRIu64",%f\n", (*it).p1, (*it).p2, (*it).p2 - (*it).p1, norm_score);
                fprintf(out_results, "%" PRIu64",%" PRIu64",%" PRIu64",%f\n", (*it).p1, (*it).p2, (*it).p2 - (*it).p1, norm_score);
            }
            
            //getchar();

            //if (i == 100) break;
            ++i;

        }
        std::cout << "\n";

        write_regions(d, query_l, selected_hash, raw_out);
        write_average_regions(d, query_l, selected_hash, average_out);

        fclose(raw_out); fclose(density_out); fclose(average_out);
        if(out_results != NULL && out_results != stdout) fclose(out_results);

    }

    // @@@@@@@@@@@@@@@@@@@@@@@ Checking hotspots

    if(hotspots != NULL){

        std::list<Hotspot> * hotspot_list = load_hotspots(hotspots);
        
        check_density_hotspots(d, frequency_table, hotspot_list, query_l, selected_hash);

        hotspot_list->sort(compare_hotspots);

        
        std::list<Hotspot>::iterator hit;
        i = 0;
        std::cout << "Start" << "\t\t" << "End" << "\t\t" << "Mean" << "\t\t" << "Ratio" << "\n";
        for (hit=hotspot_list->begin(); hit!=hotspot_list->end(); ++hit) {

            if(i<50) std::cout << (*hit).start << "\t\t" << (*hit).end << "\t\t" << (*hit).score << "\t\t" << (*hit).score/(float) theta << "\t\t" << 1/(*hit).score << "/" << 1/theta << "\n";
            ++i;

        }
        

    }


    fclose(query);
    if(hotspots != NULL) fclose(hotspots);

    free(frequency_table);
    free(d);
    free(exp_name);

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

char * init_args(int argc, char ** av, FILE ** query, FILE ** hotspots, FILE ** out_results, uint64_t * custom_kmer, uint64_t * selected_hash){

    int pNum = 0;
    char * exp_name = (char *) malloc(MAX_NAME*sizeof(char)); if(exp_name == NULL) terror("Could not generate output file names");
    uint64_t word_l = 0;
    *out_results = stdout;
    while(pNum < argc){
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           freqgen -query [query]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -out        [outfile]   (default stdout)\n");
            fprintf(stdout, "           -kmer       [Integer:   k>1 (default 2)]\n");
            fprintf(stdout, "           -word       [String:    word to plot (default AA)\n");
            fprintf(stdout, "           -lambda     [Float:     Scalar for minimum average threshold\n");
            fprintf(stdout, "           -hotspots   [File]      File containing hotspots\n");
            fprintf(stdout, "           --help      Shows help for program usage\n");
            fprintf(stdout, "\n");
            exit(1);
        }
        else if(strcmp(av[pNum], "-query") == 0){
            *query = fopen(av[pNum+1], "rt");
            if(query==NULL) terror("Could not open query file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-hotspots") == 0){
            *hotspots = fopen(av[pNum+1], "rt");
            if(hotspots==NULL) terror("Could not open hotspots file");
            pNum+=2;
        }
        else if(strcmp(av[pNum], "-out") == 0){
            exp_name[0] = '\0';
            strcpy(exp_name, av[pNum+1]);
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

    return exp_name;
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

std::list<Hotspot> * load_hotspots(FILE * hotspots){
    
    char buffer[BUFFER];
    std::list<Hotspot> * hotspot_list = new std::list<Hotspot>();

    while(!feof(hotspots)){

        if(fgets(buffer, BUFFER, hotspots) != NULL){
            
            Hotspot h;
            h.score = 0;
            sscanf(buffer,"%*s %" PRIu64" %" PRIu64, &h.start, &h.end);

            //std::cout << h.start << " " << h.end << "\n";
            //getchar();

            hotspot_list->push_back(h);

        }
    }

    return hotspot_list;
}

void generate_distribution(uint64_t * frequency_table, uint64_t * frequency_nucl, Dictionary * d, uint64_t k, FILE * query, uint64_t query_l){

    
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

            if(c == 'A') { hash = (hash << 2) + 0; ++frequency_nucl[0]; }
            if(c == 'C') { hash = (hash << 2) + 1; ++frequency_nucl[1]; }
            if(c == 'G') { hash = (hash << 2) + 2; ++frequency_nucl[2]; }
            if(c == 'T') { hash = (hash << 2) + 3; ++frequency_nucl[3]; }

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
            d[current_len - k].position = (int64_t) (current_len - k);

            // Overlapping
            //memmove(&curr_kmer[0], &curr_kmer[1], k-1);
            //--word_size;
        }
    }

    free(seq_x);

}

void generate_plotting_vector(uint64_t * frequency_table, uint64_t * frequency_nucl, Dictionary * d, uint64_t k, uint64_t query_l, uint64_t selected_hash, FILE * out_results, uint64_t factor_scale){

    // Format is:
    // desired string (e.g. ACAC)
    // frequency of the desired string (e.g. 1021321)
    // length of the sequence (e.g. 15500000)
    // probability of string given nucleotide frequencies (e.g. 0.05)
    // kmers @ pixel 1 (e.g. 10)
    // kmers @ pixel 2
    // etc.    
    // kmers @ pixel 1 but smoothed
    // kmers @ pixel 2 but smoothed
    // etc.

    char word[k+1];
    perfect_hash_to_word(word, selected_hash, k);
    word[k] = '\0';
    float p_c = 1; //observed probability
    uint64_t s;
    for(s=0; s<k; s++){
        if(word[s] == 'A') p_c = p_c * (float)frequency_nucl[0]/(float)query_l; 
        if(word[s] == 'C') p_c = p_c * (float)frequency_nucl[1]/(float)query_l;
        if(word[s] == 'G') p_c = p_c * (float)frequency_nucl[2]/(float)query_l;
        if(word[s] == 'T') p_c = p_c * (float)frequency_nucl[3]/(float)query_l;
    }
    

    fprintf(out_results, "%s\n", word);
    fprintf(out_results, "%" PRIu64"\n", frequency_table[selected_hash]);
    fprintf(out_results, "%" PRIu64"\n", query_l);
    fprintf(out_results, "%f\n", p_c);

    uint64_t * bins = (uint64_t *) calloc(VECTOR_L, sizeof(uint64_t));
    if(bins == NULL) terror("Could not allocate bins vector");

    /*
    float * bins_ow = (float *) calloc(VECTOR_L, sizeof(float));
    if(bins_ow == NULL) terror("Could not allocate bins 2 vector");

    float * fkt = (float * ) calloc(query_l, sizeof(float));
    if(fkt == NULL) terror("Could not allocate fkt vector");

    uint64_t x, y, x_sum, x_total;
    for(x=0; x<query_l; x++){
        x_sum = 0;
        x_total = 0;
        for(y=0; y<1000; y++){
            if(x+y < query_l){
                if(d[x+y].hash == selected_hash && d[x+y].position != -1){
                    ++x_sum;
                }
                ++x_total;
            }
        }
        fkt[x] = ((float)x_sum/(float)x_total);
        //printf("%" PRIu64" from x_sum : %" PRIu64"\n", fkt[x], x_sum);
    }

    for(x=0; x<VECTOR_L; x++){
        for(y=x*factor_scale; y<(x+1)*factor_scale; y++){
            if(y < query_l) bins_ow[x] += fkt[y];
        }
    }
    */
    
    uint64_t i, j = 0, current;
    for(i=0; i<VECTOR_L; i++){
        
        current = 0;

        while(d[j].position < (int64_t) ((i+1)*factor_scale) && ( d[j].position == -1 || d[j].position >= (int64_t) (i*factor_scale) ) && j < query_l){

            if(d[j].hash == selected_hash && d[j].position != -1){
                ++current;
            }

            if(current > factor_scale) terror("Higher than current");

            ++j;
            
        }
        bins[i] = current;
        //fprintf(out_results, "%" PRIu64"\n", current);

    }

    uint64_t * smoothed = average_smooth(VECTOR_L, bins, 10);
    for(i=0; i<VECTOR_L; i++) fprintf(out_results, "%" PRIu64"\n", bins[i]);
    for(i=0; i<VECTOR_L; i++) fprintf(out_results, "%" PRIu64"\n", smoothed[i]);
    //for(i=0; i<VECTOR_L; i++) fprintf(out_results, "%f\n", bins_ow[i]);

    free(bins);
    free(smoothed);
    //free(bins_ow);
    //free(fkt);
}

std::list<Region> * detect_density_regions(Dictionary * d, uint64_t * frequency_table, uint64_t query_l, uint64_t selected_hash, FILE * out_results){

    theta = (long double) frequency_table[selected_hash] / (long double) query_l;
    long double theta_inv = 1 / theta;
    fprintf(stdout, "[INFO] Theta = %Le, period = %Le\n", theta, theta_inv);
    uint64_t i = 0, j, p_i, end_position = 0;
    long double score = 0, highest_score = 0;

    std::list<Region> * density_list = new std::list<Region>();


    while(i<query_l){

        //fprintf(stdout, "%" PRIu64"\r", i);

        if(d[i].position != -1 && d[i].hash == selected_hash){

            Region r;

            p_i = 1;
            end_position = i;
            highest_score = 0;

            for(j=i+1; j<query_l; j++){

                if(d[j].position != -1 && d[j].hash == selected_hash) ++p_i;

                score = ((long double) p_i) / (long double) (j-i) - (long double)THRESHOLD * theta;

                //fprintf(stdout, "p_i: %" PRIu64" j-i: %" PRIu64" -> score: %Le\n", p_i, j-i, score);
                //getchar();

                if(score > highest_score || (j-i) < (uint64_t) theta_inv){
                    highest_score = score;
                    end_position = j;
                }

                if(score < 0.0 ) break;

            }

            if(highest_score > 0.0){

                r.p1 = i;
                r.p2 = end_position;
                r.score = (float) highest_score;

                //std::cout << "Pushing " << r.p1 << " " << r.p2 << " with p_i " << p_i << "/" << r.p2-r.p1 << " -> " << r.score << "\n";
                //getchar();

                density_list->push_back(r);
                i = end_position;

            }

            //fprintf(stdout, "@%" PRIu64" to %" PRIu64" -> S:%Le\n", i, end_position, highest_score);
            
        }

        ++i;
    }

    return density_list;
}

void check_density_hotspots(Dictionary * d, uint64_t * frequency_table, std::list<Hotspot> * hotspot_list, uint64_t query_l, uint64_t selected_hash){

    theta = (long double) frequency_table[selected_hash] / (long double) query_l;
    long double theta_inv = 1 / theta;
    fprintf(stdout, "[INFO] Theta = %Le, period = %Le\n", theta, theta_inv);
    uint64_t i = 0, p_i = 0;

    std::list<Hotspot>::iterator it;

    for (it=hotspot_list->begin(); it!=hotspot_list->end(); ++it) {

        i = (*it).start;

        while(i<(*it).end){

            p_i = 0;

            if(d[i].position != -1 && d[i].hash == selected_hash) ++p_i;

            ++i;

        }
        
        (*it).score = ((long double) p_i) / (long double) ((*it).end - (*it).start);

    }

}

// This sort is inverted @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

bool compare_regions(const Region& first, const Region& second)
{
    /*
    if(first.score == second.score) return (first.p1 > second.p1);
    return (first.score > second.score);
    */
    if(first.p2-first.p1 == second.p2-second.p1) return (first.score > second.score);
    return (first.p2-first.p1 > second.p2-second.p1);
}

bool compare_hotspots(const Hotspot& first, const Hotspot& second)
{
    /*
    if(first.score == second.score) return (first.p1 > second.p1);
    return (first.score > second.score);
    */
    if(first.score == second.score) return (first.start < second.start);
    return (first.score > second.score);
}

void write_regions(Dictionary * d, uint64_t query_l, uint64_t selected_hash, FILE * out_results){
    uint64_t i;
    for(i=0; i<query_l; i++){
        if(d[i].position != -1 && d[i].hash == selected_hash){
            fprintf(out_results, "1\n");
        }else{
            fprintf(out_results, "0\n");
        }
    }
}

void write_average_regions(Dictionary * d, uint64_t query_l, uint64_t selected_hash, FILE * out_results){
    uint64_t i, j;
    uint64_t sum, total;
    for(i=0; i<query_l; i++){
        sum = 0;
        total = 0;
        for(j=i; j<i+100; j++){
            if(d[j].position != -1 && d[j].hash == selected_hash){
                ++sum;
            }
            ++total;
        }
        fprintf(out_results, "%f\n", (float)sum/float(total));
        
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