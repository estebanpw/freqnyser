#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>

//Structs required for the dotplot workflow
#define MAXLID 200
//#define READBUF 2000000000 //2 GB
#define READBUF 500000000 //50MB
#define INITSEQS 3000 //Number of starting sequences (in database)
#define POINT 4

#define FIXED_K 12
#define TOTAL_ENTRIES 16777216

#define MAXLID 200
#define ALIGN_LEN 60 //For NW alignment
#define MAX_READ_SIZE 20000 //Maximum length of read to have a portion of the table allocated
#define MAX_WINDOW_SIZE 500 //Maximum window length to explore NW table
//#define POOL_SIZE 2500000000 // by 16 bytes it is 40 GB
#define POOL_SIZE 12500000 // 1 GB if 16 bytes
#define MAX_MEM_POOLS 256 

#define BYTES_IN_MER 4
#define MAX_DECOMP_HASH 10
#define UNSET 0
#define FALSE 1
#define TRUE 2

#define FORWARD 0
#define REVERSE 1

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))


typedef struct dictionary{
    uint64_t hash;
    int64_t position;
} Dictionary;

#endif
