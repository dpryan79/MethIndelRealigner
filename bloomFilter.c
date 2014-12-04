#include "realigner.h"
#include "murmur3.h"
#include <math.h>
int MAXLEVELS = 4;

void setMAXLEVELS(int max) {
    MAXLEVELS = max;
}

//We could also do this once with a max width
bf ** bf_init(int32_t width, int kmer) {
    uint32_t n = 10*((width+4*kmer)-1);
    //Number of bits for 10% FP rate=-n/log(0.9)
    uint64_t mask = llrint(-1.0*n/log(0.9));
//    bf **filt = (bf*) malloc(sizeof(bf*));
    bf **filt = malloc(sizeof(bf*) * MAXLEVELS);
    int i;
    assert(filt);

    for(i=0; i<MAXLEVELS; i++) {
        filt[i] = (bf*) malloc(sizeof(bf*));
        assert(filt[i]);
    }

    //Get the next highest power of 2 for the mask
    mask -= 1;
    mask |= mask>>1;
    mask |= mask>>2;
    mask |= mask>>4;
    mask |= mask>>8;
    mask |= mask>>16;
    mask |= mask>>32;
    mask += 1;
    if(mask>0x4000000000000000) mask=0x7FFFFFFFFFFFFFFF; //mask1-mask2 should fit in int64_t

    //len is the next highest multiple of 8
    filt[0]->mask = mask-1;
    filt[0]->len = (mask>>3) + ((mask&7) ? 1:0);
#ifdef DEBUG
    fprintf(stderr, "[bf_init] Creating bloom filter of size %" PRId64 " with mask %" PRIu64 "\n", filt[0]->len, filt[0]->mask);
    fflush(stderr);
#endif
    for(i=0; i<MAXLEVELS; i++) {
        filt[i]->bf = calloc(filt[0]->len, sizeof(uint8_t));
        assert(filt[i]->bf);
    }

    return filt;
}

void bf_destroy(bf **bf) {
    int i;
    for(i=0; i<MAXLEVELS; i++) {
        free(bf[i]->bf);
        free(bf[i]);
    }
    free(bf);
}

//Not used, but could be if we used a fixed filter size
void bf_reset(bf **bf) {
    int i, j;
    for(j=0; j<MAXLEVELS; j++) {
        for(i=0; i<bf[j]->len; i++) bf[j]->bf[i] = 0;
    }
}

//Add a uint64_t hash to the bloom filter
void bf_add(bf **bf, uint64_t hash) {
    uint64_t masked = hash & bf[0]->mask;
    uint8_t bit = 128>>(hash % 8ULL);
    int level = bf_exists(bf, hash);
    if(level<MAXLEVELS) bf[level]->bf[masked>>3] |= bit;
}

//Like bf_add, but this will always add the hash to the bottom-most bloom filter
void bf_add_bottom(bf **bf, uint64_t hash) {
    uint64_t masked = hash & bf[0]->mask;
    uint8_t bit = 128>>(hash % 8ULL);
    bf[MAXLEVELS-1]->bf[masked>>3] |= bit;
}

//If the bit isn't set, return 0, otherwise, return the next level with it unset (up to MAXLEVELS)
inline int bf_exists(bf **bf, uint64_t hash) {
    uint64_t masked = hash & bf[0]->mask;
    uint8_t bit = 128>>(hash % 8ULL);
    int i;
    for(i=0; i<MAXLEVELS; i++) {
        if(bf[i]->bf[masked>>3] & bit) continue;
        return i;
    }
    return i;
}

//Like bf_exists, but only checks the bottom level
inline int bf_exists_bottom(bf **bf, uint64_t hash) {
    uint64_t masked = hash & bf[0]->mask;
    uint8_t bit = 128>>(hash % 8ULL);
    return (bf[MAXLEVELS-1]->bf[masked>>3] & bit);
}

//We only support up to a 64bit hash, even though 128 is computed
uint64_t hash_seq(char *seq, int len) {
    uint64_t hash_val[2];
    uint32_t seed = 0xAAAAAAAA;
#if UINTPTR_MAX == 0xffffffff
    MurmurHash3_x86_128((void *) seq, len, seed, (void *) &hash_val);
#else
    //Assuming 64bit is probably not terrible these days
    MurmurHash3_x64_128((void *) seq, len, seed, (void *) &hash_val);
#endif
    return hash_val[0];
}
