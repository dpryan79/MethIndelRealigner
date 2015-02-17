#include "realigner.h"
#include "murmur3.h"
#include <string.h>
#include <assert.h>

/*******************************************************************************
*
*  This currently uses an unordered linked-list, which may not perform optimally
*
*******************************************************************************/

hashTable * ht_init(int32_t nAlignments, int threshold) {
    uint64_t n = 500+2*((uint64_t) nAlignments);
    hashTable *ht = malloc(sizeof(hashTable));
    assert(ht);

    ht->threshold = threshold;
    ht->n = n;
    ht->entries = calloc(n, sizeof(hashTableEntry *));
    assert(ht->entries);
    return ht;
}

hashTableEntry * ht_makeEntry(char *seq, int l) {
    hashTableEntry *e = calloc(1, sizeof(hashTableEntry));
    assert(e);

    e->cnt = 1;
    e->seq = strndup(seq, l);
    assert(e->seq);

    return e;
}

hashTableEntry * ht_hasEntry(hashTable *ht, uint64_t h, char *seq, int len) {
    uint64_t n = h%ht->n;
    hashTableEntry *p = ht->entries[n];

    while(p) {
        if(strncmp(p->seq,seq, len) == 0) break;
        p = p->next;
    }
    return p;
}

void ht_addEntry(hashTable *ht, uint64_t h, hashTableEntry *e) {
    uint64_t n = h%ht->n;
    hashTableEntry *p = ht->entries[n];

    if(!p) {
        ht->entries[n] = e;
    } else {
        while(p->next) p=p->next;
        p->next = e;
    }
}

void ht_addIncrement(hashTable *ht, char *seq, int len) {
    uint64_t h = hash_seq(seq, len);
    hashTableEntry *e;
    e = ht_hasEntry(ht, h, seq, len);

    if(e) {
        if(e->cnt < 0xFF) {
            e->cnt++;
        }
    } else {
        e = ht_makeEntry(seq, len);
        ht_addEntry(ht, h, e);
    }
}

void ht_addIncrementMax(hashTable *ht, char *seq, int len) {
    uint64_t h = hash_seq(seq, len);
    hashTableEntry *e;
    e = ht_hasEntry(ht, h, seq, len);

    if(e) e->cnt = 0xFF;
    else {
        e = ht_makeEntry(seq, len);
        e->cnt = 0xFF;
        ht_addEntry(ht, h, e);
    }
}

int ht_sufficientCount(hashTable *ht, char *seq, int len) {
    uint64_t h = hash_seq(seq, len);
    hashTableEntry *e;
    e = ht_hasEntry(ht, h, seq, len);

    if(!e) return 0;
    if(e->cnt < ht->threshold) return 0;
    return 1;
}

uint64_t ht_numEntries(hashTable *ht) {
    uint64_t i, n = 0;
    hashTableEntry *p;

    for(i=0; i<ht->n; i++) {
        p = ht->entries[i];
        while(p) {
            n++;
            p = p->next;
        }
    }
    return n;
}

uint64_t ht_maxDepth(hashTable *ht) {
    uint64_t i, n, biggestN = 0;
    hashTableEntry *p;

    for(i=0; i<ht->n; i++) {
        n = 0;
        p = ht->entries[i];
        while(p) {
            n++;
            p = p->next;
        }
        if(n>biggestN) biggestN = n;
    }
    return biggestN;
}

void ht_destroy(hashTable *ht) {
    uint64_t i;
    hashTableEntry *prev, *next;

    for(i=0; i<ht->n; i++) {
        prev = ht->entries[i];
        while(prev) {
            next = prev->next;
            free(prev->seq);
            free(prev);
            prev = next;
        }
    }
    free(ht->entries);
    free(ht);
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

