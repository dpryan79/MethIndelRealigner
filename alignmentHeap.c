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
alignmentHeap * writeHeapUntil(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap) {
    alignmentHeap *new = alignmentHeap_init(1000); //CHANGE ME!
    int cur_i = 0;
    bam1_t *b = heap->heap[0];
    int32_t lpos = -1;
    InDel reg;

    if(lastTargetNode == NULL) {
        writeHeap(of, hdr, heap);
        return new;
    }
    if(lastTargetNode->tid != heap->heap[0]->core.tid) {
        writeHeap(of, hdr, heap);
        return new;
    }
    reg.tid = b->core.tid;
    reg.start = b->core.pos;
    reg.end = bam_endpos(b)+1;
    while(TargetNodeCmp(&reg, lastTargetNode) != 0) {
        sam_write1(of, hdr, b);
        bam_destroy1(heap->heap[cur_i++]);
        //We hit the end of the heap before the next ROI
        if(cur_i >= heap->l) {
            free(heap->heap);
            free(heap);
            return new;
        }
        b = heap->heap[cur_i];
        reg.tid = b->core.tid;
        reg.start = b->core.pos;
        reg.end = bam_endpos(b)+1;
    }
    lpos = b->core.pos;
    //Keep sort order, processing reads with same start pos
    while(b->core.pos == lpos) {
        if(TargetNodeCmp(&reg, lastTargetNode) == 0) {
            new->heap[new->l++] = b;
        } else {
            sam_write1(of, hdr, b);
            bam_destroy1(b);
        }
        if(++cur_i >= heap->l) {
            free(heap->heap);
            free(heap);
            return new;
        }
        b = heap->heap[cur_i];
        reg.tid = b->core.tid;
        reg.start = b->core.pos;
        reg.end = bam_endpos(b)+1;
    }
    //Transfer the remainder over
    for(;cur_i < heap->l; cur_i++) new->heap[new->l++] = heap->heap[cur_i];
    free(heap->heap);
    free(heap);

    return new;
}
