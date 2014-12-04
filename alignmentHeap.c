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
void writeHeap(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, int doSort) {
    int i;
    if(doSort) sortHeap(heap);
    for(i=0; i<heap->l; i++) {
        sam_write1(of, hdr, heap->heap[i]);
        bam_destroy1(heap->heap[i]);
    }
}

alignmentHeap * mergeHeap(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, bam1_t *curb, samFile *fp) {
    int cur_i = 0;

    while(cur_i < heap->l && (heap->heap[heap->l-1]->core.pos > curb->core.pos && heap->heap[heap->l-1]->core.tid == curb->core.tid)) {
        if(heap->heap[cur_i]->core.pos <= curb->core.pos) {
            sam_write1(of, hdr, curb);
            bam_destroy1(heap->heap[cur_i++]);
        } else {
            sam_write1(of, hdr, curb);
            if(sam_read1(fp, hdr, curb) <= 0) break; //EOF
        }
    }
    //In case we hit EOF before the end of the heap...
    while(cur_i < heap->l) {
        sam_write1(of, hdr, heap->heap[cur_i]);
        bam_destroy1(heap->heap[cur_i++]);
    }
    heap->l = 0;
    return heap;
}

//We could get a bit more elaborate and 
alignmentHeap * mergeHeapROI(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, bam1_t *curb, samFile *fp) {
    int cur_i = 0, status = 0, newl = 0;

    //while(cur_i < heap->l && (heap->heap[heap->l-1]->core.pos > curb->core.pos && heap->heap[heap->l-1]->core.tid == curb->core.tid)) {
    while(cur_i < heap->l) {
        if(heap->heap[cur_i]->core.pos <= curb->core.pos) {
            if(bam_endpos(heap->heap[cur_i])-1 < lastTargetNode->start || heap->heap[cur_i]->core.tid != lastTargetNode->tid) {
                sam_write1(of, hdr, heap->heap[cur_i]);
                bam_destroy1(heap->heap[cur_i++]);
            } else {
                status = 1;
                break;
            }
        } else {
            if(bam_endpos(curb)-1 < lastTargetNode->start && curb->core.tid <= lastTargetNode->tid) {
                sam_write1(of, hdr, curb);
                if(sam_read1(fp, hdr, curb) <= 0) { //EOF
                    status = 2;
                    break;
                }
            } else {
                status = 1;
                break;
            }
        }
    }

    if(status == 0) { //Merged region before next ROI, 
        assert(cur_i == heap->l);
        heap->l = 0;
    } else if(status == 1) { //Something is overlapping!
        if(cur_i < heap->l) { //Shift everything
            for(;cur_i < heap->l; cur_i++) heap->heap[newl++] = heap->heap[cur_i];
            heap->l = newl;
        } else {
            heap->l = 0;
        }
    } else { //Hit EOF, the remaining heap should be processed
        //This shouldn't actually happen unless someone uses a target file for a different BAM!
        assert(1==0);
    }
    return heap;
}

//Export
alignmentHeap * writeHeapUntil(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, bam1_t *curb, samFile *fp) {
    int cur_i = 0;
    bam1_t *b;
    int32_t lpos = -1, newl = 0;
    InDel reg;

    sortHeap(heap);

    /****************************************
    * If we moved a start position to the right of an ROI, then the output may become mis-sorted.
    * We end up running into two cases:
    * (1) The one of the alignments in a mis-sorted region overlaps the next ROI
    * (2) The mis-sorted region can be merged before hitting the next ROI.
    * Option 2 is easy enough to deal with, while option 1 is more annoying.
    ****************************************/
    if(curb->core.tid == heap->heap[0]->core.tid && curb->core.pos < heap->heap[heap->l-1]->core.pos) { //Need to merge!
        if(lastTargetNode == NULL) { //No ROI to care about, merge until finished
            heap = mergeHeap(of, hdr, heap, curb, fp);
        } else { //Merge until we hit an ROI
            heap = mergeHeapROI(of, hdr, heap, curb, fp);
        }
        return heap;
    }

    if(lastTargetNode == NULL) {
        writeHeap(of, hdr, heap, 0);
        heap->l = 0;
        return heap;
    }
    heap->start = lastTargetNode->start;
    heap->end = lastTargetNode->end;
    if(lastTargetNode->tid != heap->heap[0]->core.tid) {
        writeHeap(of, hdr, heap, 0);
        heap->l = 0;
        return heap;
    }

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
