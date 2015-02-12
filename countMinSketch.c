#include "realigner.h"
#include "murmur3.h"
#include <math.h>
#define MAXCNT 15

//Create an empty count-min sketch
//width is the width of the realigned region
cms * cms_init(int32_t width, int kmer, int threshold) {
    uint32_t n = 5*((width+4*kmer)-1); //In an ideal world, we'd use hyperloglog...
    //Number of entries for 10% FP rate=-n/log(0.9)
    uint64_t mask = llrint(-1.0*n/log(0.9));
    cms *cmsOut = calloc(1, sizeof(cms));
    assert(cmsOut);

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
    cmsOut->threshold = threshold;
    cmsOut->mask = mask-1;
    cmsOut->len = (mask>>1) + 1; //4 bit counter per entry, so 1 byte per 2 entries
    cmsOut->cnts = calloc(cmsOut->len, sizeof(uint8_t));
    assert(cmsOut->cnts);

#ifdef DEBUG
    fprintf(stderr, "[cms_init] Creating a count-min sketch of size %" PRId64 " with mask %" PRIu64 " and threshold %i\n", cmsOut->len, cmsOut->mask, cmsOut->threshold);
    fflush(stderr);
#endif

    return cmsOut;
}

void cms_destroy(cms *cntMS) {
    free(cntMS->cnts);
    free(cntMS);
}

//Increment a count in the count-min sketch
void cms_increment(cms *cntMS, uint64_t hash) {
    uint64_t masked = hash & cntMS->mask;
    uint8_t lower = cntMS->cnts[masked>>1] & 0xF;
    uint8_t upper = (cntMS->cnts[masked>>1] & 0xF0)>>4;
#ifdef DEBUG
    fprintf(stderr, "[cms_increment] Incrementing count of %"PRIu64" from %"PRIu8" to ", hash, (masked&1)?lower:upper);
#endif
    if(masked&1) {
        if(lower<MAXCNT) lower++;
    } else {
        if(upper<MAXCNT) upper++;
    }
#ifdef DEBUG
    fprintf(stderr, "%"PRIu8"\n", (masked&1)?lower:upper);
#endif
    cntMS->cnts[masked>>1] = (upper<<4) | lower;
}

//Like cms_increment, but will always max a value
void cms_max(cms *cntMS, uint64_t hash) {
    uint64_t masked = hash & cntMS->mask;
    uint8_t lower = cntMS->cnts[masked>>1] & 0xF;
    uint8_t upper = (cntMS->cnts[masked>>1] & 0xF0)>>4;
    if(masked&1) lower = MAXCNT;
    else upper = MAXCNT;
    cntMS->cnts[masked>>1] = (upper<<4) | lower;
}

//Return the counter for a given hash
inline int cms_val(cms *cntMS, uint64_t hash) {
    uint64_t masked = hash & cntMS->mask;
    if(masked&1) return cntMS->cnts[masked>>1] & 0xF;
    return (cntMS->cnts[masked>>1] & 0xF0)>>4;
}

//Like cms_val, but checks that the cms is at least the minimum kmer count
inline int cms_val_sufficient(cms *cntMS, uint64_t hash) {
    return (cms_val(cntMS, hash) >= cntMS->threshold);
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
