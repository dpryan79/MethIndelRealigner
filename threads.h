#include <semaphore.h>
#ifndef REALIGN_H
#include "realigner.h"
#endif

//alignReads2Paths
pthread_t *alignReads2Paths_threads;
sem_t alignReads2Paths_sem1, alignReads2Paths_sem2;
pthread_spinlock_t alignReads2Paths_spinlock, alignReads2Paths_finished_spinlock;

//SemiGlobalAlignment
pthread_t *SemiGlobalAlignment_threads;
sem_t SemiGlobalAlignment_sem1, SemiGlobalAlignment_sem2;
pthread_spinlock_t SemiGlobalAlignment_spinlock, SemiGlobalAlignment_finished_spinlock;

alignmentHeap *GLOBAL_heap;
int GLOBAL_k, GLOBAL_nt, GLOBAL_len, GLOBAL_finished, *thread_ids;
int32_t GLOBAL_start, GLOBAL_end;
bf *GLOBAL_filter;
char *GLOBAL_CT, *GLOBAL_GA, *GLOBAL_refSeq;
int8_t *GLOBAL_ref;
int32_t GLOBAL_refLen, GLOBAL_maxIns, GLOBAL_maxDel, GLOBAL_refIndex;
int32_t *GLOBAL_readLBound, *GLOBAL_readRBound, GLOBAL_refLBound, GLOBAL_refRBound;
s_align ***GLOBAL_salign, **GLOBAL_salign2;
paths *GLOBAL_CTpaths, *GLOBAL_GApaths, *GLOBAL_paths;

s_align ***signal_alignReads2Paths(int32_t *readLBound, int32_t *readRBound);
s_align **signal_SemiGlobalAlignment(int32_t *refIndex);
void *threaded_alignReads2Paths(void *a);
void *threaded_SemiGlobalAlignment(void *a);
void setup_alignReads2Paths(s_align ***readal, alignmentHeap *heap, paths *CTpaths, paths *GApaths, int32_t refLBound, int32_t refRBound, int32_t *readLBound, int32_t *readRBound, int k);
void setup_SemiGlobalAlignment(s_align **alignments, int8_t *ref, int32_t refLen, paths *p, int32_t k, char *refSeq);
void setup_threads(int nt);
void destroy_threads();
