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

int alignmentHeapSortFunc(const void *pa, const void *pb) {
    bam1_t *a = *((bam1_t **) pa);
    bam1_t *b = *((bam1_t **) pb);
    int rv = 0;

    if(a->core.tid == b->core.tid) {
        if(a->core.pos == b->core.pos) {
            if(bam_endpos(a) != bam_endpos(b)) rv = bam_endpos(a) - bam_endpos(b);
        } else {
            rv = a->core.pos - b->core.pos;
        }
    } else {
        rv = a->core.tid - b->core.tid;
    }
    return rv;
}

void sortHeap(alignmentHeap *heap) {
    qsort((void *) heap->heap, heap->l, sizeof(bam1_t*), &alignmentHeapSortFunc);
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
    sortHeap(heap);
    for(i=0; i<heap->l; i++) {
        sam_write1(of, hdr, heap->heap[i]);
        bam_destroy1(heap->heap[i]);
    }
}

//Export
alignmentHeap * writeHeapUntil(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, int depth) {
    int cur_i = 0;
    bam1_t *b;
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

    sortHeap(heap);
    b = heap->heap[0];
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
    }
    //Transfer the remainder over
    for(;cur_i < heap->l; cur_i++) {
        heap->heap[newl++] = heap->heap[cur_i];
    }
    heap->l = newl;
#ifdef DEBUG
    fprintf(stderr, "[writeHeapUntil] Returning a heap of length %" PRId32 "\n", newl); fflush(stderr);
#endif

    return heap;
}
