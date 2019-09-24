#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H
#include "structs.h"
/**
 * Print the error message 's' and exit(-1)
 */
void terror(const char *s);


/**
 * Function to read char by char buffered from a FILE
 */
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

uint64_t quick_pow4(uint64_t n);

uint64_t quick_pow4byLetter(uint64_t n, const char c);

uint64_t hashOfWord(const unsigned char * word, uint32_t k, uint64_t offset);

void perfect_hash_to_word(char * word, uint64_t hash, uint32_t k);

uint64_t collisioned_hash(const unsigned char * word, uint32_t k);

void decomposed_hash_of_word(const unsigned char * word, unsigned char * vector, uint32_t k);

uint64_t xor_decomposed_hash(unsigned char * vector1, unsigned char * vector2, uint32_t k);

uint64_t asciiToUint64(const char *text);

unsigned char complement(unsigned char c);

void inplace_reverse_and_complement(unsigned char *d, uint64_t l);

uint64_t * average_smooth(uint64_t l, uint64_t * v, uint64_t window);


#endif /* COMMON_FUNCTIONS_H */
