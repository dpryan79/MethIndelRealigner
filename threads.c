#include <inttypes.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include "threads.h"

//external functions
extern s_align **alignReads2Paths(bam1_t *b, int strand, int32_t *subreadM, int8_t **subreadSeq, paths *p, int32_t refLBound, int32_t refRBound, int32_t *readLBound, int32_t *readRBound, int k);
extern int getStrand(bam1_t *b);
extern s_align **alignPaths2Ref(paths *p, int8_t *ref, int32_t refLen, char *refSeq, int32_t *refIndex, int32_t k);

s_align ***signal_alignReads2Paths(int32_t *readLBound, int32_t *readRBound) {
    int i;

    pthread_spin_lock(&alignReads2Paths_finished_spinlock); //This prevents a race
    pthread_spin_unlock(&alignReads2Paths_spinlock); //Release the workers

    //Wait until they're finished
    for(i=0; i<GLOBAL_nt; i++) sem_wait(&alignReads2Paths_sem1);

    //Release the loop mutex, ensure the workers can loop before returning
    pthread_spin_lock(&alignReads2Paths_spinlock); //prevent a worker from restarting
    pthread_spin_unlock(&alignReads2Paths_finished_spinlock); //Release the workers to loop
    for(i=0; i<GLOBAL_nt; i++) sem_wait(&alignReads2Paths_sem2); //Don't race the workers to the finished spinlock

    readLBound = GLOBAL_readLBound;
    readRBound = GLOBAL_readRBound;
    return GLOBAL_salign;
}

s_align ** signal_SemiGlobalAlignment(int32_t *refIndex) {
    int i;

    pthread_spin_lock(&SemiGlobalAlignment_finished_spinlock); //This prevents a race
    pthread_spin_unlock(&SemiGlobalAlignment_spinlock); //Release the workers

    //Wait until they're finished
    for(i=0; i<GLOBAL_nt; i++) sem_wait(&SemiGlobalAlignment_sem1);

    //Release the loop mutex, ensure the workers can loop before returning
    pthread_spin_lock(&SemiGlobalAlignment_spinlock); //prevent a worker from restarting
    pthread_spin_unlock(&SemiGlobalAlignment_finished_spinlock); //Release the workers to loop
    for(i=0; i<GLOBAL_nt; i++) sem_wait(&SemiGlobalAlignment_sem2); //Don't race the workers to the finished spinlock

    *refIndex = GLOBAL_refIndex;
    return GLOBAL_salign2;
}

void * threaded_alignReads2Paths(void *a) {
    int thread_id = *((int*) a);
    int i, strand;
    int8_t *subreadSeq = NULL;
    int32_t subreadM;

    subreadSeq = malloc(50*sizeof(int8_t));
    assert(subreadSeq);
    subreadM = 50;

loop:
    pthread_spin_lock(&alignReads2Paths_spinlock);
    pthread_spin_unlock(&alignReads2Paths_spinlock);

    if(GLOBAL_finished) goto cleanup;

    //The guts of the function
    for(i=thread_id; i<GLOBAL_heap->l; i+=GLOBAL_nt) {
        if(GLOBAL_heap->heap[i]->core.qual < MINMAPQ) continue;
        strand = getStrand(GLOBAL_heap->heap[i]);
        GLOBAL_salign[i] = alignReads2Paths(GLOBAL_heap->heap[i], strand, &subreadM, &subreadSeq, (strand==1) ? GLOBAL_CTpaths : GLOBAL_GApaths, GLOBAL_refLBound, GLOBAL_refRBound, GLOBAL_readLBound+i, GLOBAL_readRBound+i, GLOBAL_k);
    }

    //Lock and update
    sem_post(&alignReads2Paths_sem1); //Signal finished
    pthread_spin_lock(&alignReads2Paths_finished_spinlock); //The signalling thread acknowledges
    pthread_spin_unlock(&alignReads2Paths_finished_spinlock);
    sem_post(&alignReads2Paths_sem2); //Signal looping (otherwise, we race to lock the finished lock

    goto loop;

cleanup:
    free(subreadSeq);
    return NULL;
}
void *threaded_SemiGlobalAlignment(void *a) {
    int thread_id = *((int*) a);
    int i;

loop:
    pthread_spin_lock(&SemiGlobalAlignment_spinlock);
    pthread_spin_unlock(&SemiGlobalAlignment_spinlock);

    if(GLOBAL_finished) goto cleanup;

    //The guts of the function
    for(i=thread_id; i<GLOBAL_paths->l; i+=GLOBAL_nt) {
        GLOBAL_salign2[i] = SemiGlobalAlignment(GLOBAL_ref, GLOBAL_refLen, GLOBAL_paths->conv[i], GLOBAL_paths->len[i], GLOBAL_k);
        if(GLOBAL_refIndex == -1 && GLOBAL_salign2[i]->cigarLen == 1)
            if(strcmp(GLOBAL_refSeq, GLOBAL_paths->path[i]) == 0) GLOBAL_refIndex = i;
    }

    //Lock and update
    sem_post(&SemiGlobalAlignment_sem1); //Signal finished
    pthread_spin_lock(&SemiGlobalAlignment_finished_spinlock); //The signalling thread acknowledges
    pthread_spin_unlock(&SemiGlobalAlignment_finished_spinlock);
    sem_post(&SemiGlobalAlignment_sem2); //Signal looping (otherwise, we race to lock the finished lock

    goto loop;

cleanup:
    return NULL;
}

void setup_alignReads2Paths(s_align ***readal, alignmentHeap *heap, paths *CTpaths, paths *GApaths, int32_t refLBound, int32_t refRBound, int32_t *readLBound, int32_t *readRBound, int k) {
    GLOBAL_salign = readal;
    GLOBAL_heap = heap;
    GLOBAL_CTpaths = CTpaths;
    GLOBAL_GApaths = GApaths;
    GLOBAL_refLBound = refLBound;
    GLOBAL_refRBound = refRBound;
    GLOBAL_readLBound = readLBound;
    GLOBAL_readRBound = readRBound;
    GLOBAL_k = k;
}
void setup_SemiGlobalAlignment(s_align **alignments, int8_t *ref, int32_t refLen, paths *p, int32_t k, char *refSeq) {
    GLOBAL_salign2 = alignments;
    GLOBAL_ref = ref;
    GLOBAL_refLen = refLen;
    GLOBAL_paths = p;
    GLOBAL_k = k;
    GLOBAL_refIndex = -1;
    GLOBAL_refSeq = refSeq;
}

void setup_threads(int nt) {
    int i;
    thread_ids = malloc(sizeof(int) * nt);
    assert(thread_ids);

    for(i=0; i<nt; i++) thread_ids[i] = i;
    GLOBAL_nt = nt;
    GLOBAL_finished = 0;

    //So valgrind doesn't complain
    GLOBAL_salign = NULL;

    //alignReads2Paths
    pthread_spin_init(&alignReads2Paths_spinlock, PTHREAD_PROCESS_PRIVATE);
    pthread_spin_lock(&alignReads2Paths_spinlock);
    pthread_spin_init(&alignReads2Paths_finished_spinlock, PTHREAD_PROCESS_PRIVATE);
    sem_init(&alignReads2Paths_sem1, 0, 0);
    sem_init(&alignReads2Paths_sem2, 0, 0);
    alignReads2Paths_threads = calloc(nt, sizeof(pthread_t));
    assert(alignReads2Paths_threads);
    for(i=0; i<nt; i++) pthread_create(alignReads2Paths_threads+i, NULL, &threaded_alignReads2Paths, thread_ids+i);

    //SemiGlobalAlignment
    pthread_spin_init(&SemiGlobalAlignment_spinlock, PTHREAD_PROCESS_PRIVATE);
    pthread_spin_lock(&SemiGlobalAlignment_spinlock);
    pthread_spin_init(&SemiGlobalAlignment_finished_spinlock, PTHREAD_PROCESS_PRIVATE);
    sem_init(&SemiGlobalAlignment_sem1, 0, 0);
    sem_init(&SemiGlobalAlignment_sem2, 0, 0);
    SemiGlobalAlignment_threads = calloc(nt, sizeof(pthread_t));
    assert(SemiGlobalAlignment_threads);
    for(i=0; i<nt; i++) pthread_create(SemiGlobalAlignment_threads+i, NULL, &threaded_SemiGlobalAlignment, thread_ids+i);
}

void destroy_threads() {
    int i;
    GLOBAL_finished = 1;

    //alignReads2Paths
    pthread_spin_unlock(&alignReads2Paths_spinlock);
    for(i=0; i<GLOBAL_nt; i++) pthread_join(alignReads2Paths_threads[i], NULL);
    pthread_spin_destroy(&alignReads2Paths_spinlock);
    pthread_spin_destroy(&alignReads2Paths_finished_spinlock);
    sem_destroy(&alignReads2Paths_sem1);
    sem_destroy(&alignReads2Paths_sem2);
    free(alignReads2Paths_threads);

    //SemiGlobalAlignment
    pthread_spin_unlock(&SemiGlobalAlignment_spinlock);
    for(i=0; i<GLOBAL_nt; i++) pthread_join(SemiGlobalAlignment_threads[i], NULL);
    pthread_spin_destroy(&SemiGlobalAlignment_spinlock);
    pthread_spin_destroy(&SemiGlobalAlignment_finished_spinlock);
    sem_destroy(&SemiGlobalAlignment_sem1);
    sem_destroy(&SemiGlobalAlignment_sem2);
    free(SemiGlobalAlignment_threads);

    //Clean up
    free(thread_ids);
}
