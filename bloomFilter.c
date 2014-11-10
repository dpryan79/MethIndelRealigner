#include "realigner.h"
#include "murmur3.h"
#include <math.h>

//We could also do this once with a max width
bf * bf_init(int32_t width, int kmer) {
    uint32_t n = 4*((width+4*kmer)-1);
    //Number of bits for 10% FP rate=-n/log(0.9)
    uint64_t mask = llrint(-1.0*n/log(0.9));
    bf *filt = (bf*) malloc(sizeof(bf));
    assert(filt);

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
    filt->mask = mask-1;
    filt->len = (mask>>3) + ((mask&7) ? 1:0);
#ifdef DEBUG
    fprintf(stderr, "[bf_init] Creating bloom filter of size %" PRId64 " with mask %" PRIu64 "\n", filt->len, filt->mask);
    fflush(stderr);
#endif
    filt->bf = calloc(filt->len, sizeof(uint8_t));
    assert(filt->bf);

    return filt;
}

void bf_destroy(bf *bf) {
    free(bf->bf);
    free(bf);
}

//Not used, but could be if we used a fixed filter size
void bf_reset(bf *bf) {
    int i;
    for(i=0; i<bf->len; i++) bf->bf[i] = 0;
}

//Add a uint64_t hash to the bloom filter
void bf_add(bf *bf, uint64_t hash) {
    uint64_t masked = hash & bf->mask;
    uint8_t bit = 128>>(hash % 8ULL);
    bf->bf[masked>>3] |= bit;
}

//Return >0 if the filter may hold the hash, otherwise 0
inline int bf_exists(bf *bf, uint64_t hash) {
    uint64_t masked = hash & bf->mask;
    uint8_t bit = 128>>(hash % 8ULL);
    return (bf->bf[masked>>3] & bit);
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
