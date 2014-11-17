#include "realigner.h"
//Export
//Returns NULL on error
alignmentHeap * alignmentHeap_init(int size) {
    alignmentHeap *ret = calloc(1, sizeof(alignmentHeap));
    ret->m = size; //User configureable
    ret->heap = malloc(sizeof(bam1_t*) * size);
    if(ret->heap == NULL) {
        free(ret);
        return NULL;
    }
    return ret;
}

//Export
void alignmentHeap_destroy(alignmentHeap *heap) {
    int i;
    for(i=0; i<heap->l; i++) bam_destroy1(heap->heap[i]);
    free(heap->heap);
    free(heap);
}

//Export
void writeHeap(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap) {
    int i;
    for(i=0; i<heap->l; i++) {
        sam_write1(of, hdr, heap->heap[i]);
        bam_destroy1(heap->heap[i]);
    }
}

//Export
alignmentHeap * writeHeapUntil(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, int depth) {
    int cur_i = 0;
    bam1_t *b = heap->heap[0];
    int32_t lpos = -1, newl = 0;
    InDel reg;

    if(lastTargetNode == NULL) {
        writeHeap(of, hdr, heap);
        heap->l = 0;
        return heap;
    }
    heap->start = lastTargetNode->start;
    heap->end = lastTargetNode->end;
    if(lastTargetNode->tid != heap->heap[0]->core.tid) {
        writeHeap(of, hdr, heap);
        heap->l = 0;
        return heap;
    }
    reg.tid = b->core.tid;
    reg.start = b->core.pos;
    reg.end = bam_endpos(b)-1;
    while(reg.end < lastTargetNode->start) {
        sam_write1(of, hdr, b);
        bam_destroy1(heap->heap[cur_i++]);
        //We hit the end of the heap before the next ROI
        if(cur_i >= heap->l) {
            heap->l = 0;
            return heap;
        }
        b = heap->heap[cur_i];
        reg.tid = b->core.tid;
        reg.start = b->core.pos;
        reg.end = bam_endpos(b)-1;
    }
    lpos = b->core.pos;
    //Keep sort order, processing reads with same start pos
/*
fprintf(stderr, "[writeHeapUntil] read from %"PRId32"-%"PRId32" next starts %"PRId32" %"PRId32"\n", reg.start, reg.end, lastTargetNode->start, newl);
int i;
if(reg.end-reg.start>10000) {
    fprintf(stderr, "[writeHeapUntil] %s just got weird!\n", bam_get_qname(b));
    for(i=0; i<b->core.n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(bam_get_cigar(b)[i]), BAM_CIGAR_STR[bam_cigar_op(bam_get_cigar(b)[i])]);
    fprintf(stderr, "\n"); fflush(stderr);
}
*/
    while(b->core.pos == lpos) {
        if(reg.end >= lastTargetNode->start && reg.start < lastTargetNode->end && reg.tid == lastTargetNode->tid) {
            heap->heap[newl++] = b;
        } else {
            sam_write1(of, hdr, b);
            bam_destroy1(b);
        }
        if(++cur_i >= heap->l) { //Hit the end of the heap
            heap->l = newl;
            return heap;
        }
        b = heap->heap[cur_i];
        reg.tid = b->core.tid;
        reg.start = b->core.pos;
        reg.end = bam_endpos(b)+1;
/*
fprintf(stderr, "[writeHeapUntil] read from %"PRId32"-%"PRId32" next starts %"PRId32" %"PRId32"\n", reg.start, reg.end, lastTargetNode->start, newl);
if(reg.end-reg.start>10000) {
    fprintf(stderr, "[writeHeapUntil] %s just got weird!\n", bam_get_qname(b));
    for(i=0; i<b->core.n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(bam_get_cigar(b)[i]), BAM_CIGAR_STR[bam_cigar_op(bam_get_cigar(b)[i])]);
    fprintf(stderr, "\n"); fflush(stderr);
}
*/
    }
    //Transfer the remainder over
    for(;cur_i < heap->l; cur_i++) heap->heap[newl++] = heap->heap[cur_i];
    heap->l = newl;
#ifdef DEBUG
    fprintf(stderr, "[writeHeapUntil] Returning a heap of length %" PRId32 "\n", newl); fflush(stderr);
#endif

    return heap;
}
