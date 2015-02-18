#include "realigner.h"
#include <math.h>

extern int32_t *bam2positions(bam1_t*);
extern void findPositions(int32_t *, int, int32_t, int32_t, int *, int *);
extern char int2base[];
extern int getStrand(bam1_t *b);

kmvSketch * kmv_init(int l) {
    int i;
    kmvSketch *kmv = malloc(sizeof(kmvSketch));
    assert(kmv);

    kmv->vals = malloc(sizeof(uint64_t)*l);
    assert(kmv->vals);

    kmv->l = l;
    for(i=0; i<l; i++) {
        kmv->vals[i] = -1;
    }
    return kmv;
}

void kmv_add(kmvSketch *kmv, uint64_t v) {
    int i = kmv->l-2;
    if(kmv->vals[kmv->l-1] < v) return;

    while(kmv->vals[i] > v) {
        kmv->vals[i+1] = kmv->vals[i];
        if(i==0) break;
        i--;
    }
    kmv->vals[i] = v;
}

void kmv_addSeq(kmvSketch *kmv, char *seq, int l) {
    uint64_t h = hash_seq(seq, l);
    kmv_add(kmv, h);
}

int kmv_estimate(kmvSketch *kmv) {
    return (int) ceil(((double) (kmv->l-1))/(((double) kmv->vals[kmv->l-1])/((double) 0xFFFFFFFFFFFFFFFF)));
}

void kmv_destroy(kmvSketch *kmv) {
    free(kmv->vals);
    free(kmv);
}

//Create an array of some defined length of doubles of 1.0
//Iterate over the reference sequence k-mers
//Iterate over every Nth read, processing its k-mers
//Estimate the number of distinct values (k-1)/kmax (so (10-1)/0.1 for an array of 10 with largest value 0.1)

void bam2kmerSketch(bam1_t *b, kmvSketch *kmv, int k, int32_t start, int32_t end, char *CT, char *GA, int32_t refLen) {
    int i, start2, end2;
    int offset = (getStrand(b) & 1) ? 0 : 16;
    int32_t *positions = NULL; //Could reuse this like a kstring_t
    int32_t refStart = (start-k>=0) ? start-k : 0;
    int32_t refEnd = refStart+refLen-1;
    kstring_t *ks = calloc(1, sizeof(kstring_t));
    assert(ks);

    positions = bam2positions(b);
    findPositions(positions, b->core.l_qseq, start, end-1, &start2, &end2);
    if(start2==-1) {
        free(positions);
        return;
    }

    //Does this even overlap the ROI?
    if(positions[start2]-refStart+end2-start2+1 <= 2*k) {
        free(positions);
        return;
    }

    //Grow ks as needed
    ks_resize(ks, (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart)+1);
    ks->s[0] = '\0'; //So valgrind doesn't complain
    ks->l = (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart)+1;

    //Add 5' reference sequence
    if(positions[start2]-refStart) {
        if(offset==0) strncpy(ks->s, CT, positions[start2]-refStart);
        else strncpy(ks->s, GA, positions[start2]-refStart);
    }

    //Extract the C->T and G->A sequences
    for(i=start2; i<=end2; i++) {
        ks->s[positions[start2]-refStart+i-start2] = int2base[offset+bam_seqi(bam_get_seq(b), i)];
    }

    //Add 3' reference sequence
    if(refEnd-positions[end2]) {
        if(offset==0) strncpy(ks->s+positions[start2]-refStart+end2-start2+1,
                              CT+refLen-refEnd+positions[end2],
                              refEnd-positions[end2]);
        else strncpy(ks->s+positions[start2]-refStart+end2-start2+1,
                     GA+refLen-refEnd+positions[end2],
                     refEnd-positions[end2]);
    }
    //Add the kmers, ignoring the first and last
    for(i=1; i<ks->l-k-1; i++) kmv_addSeq(kmv, ks->s+i, k);

    free(positions);
    free(ks->s);
    free(ks);
}
