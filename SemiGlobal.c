#include "realigner.h"

void printMatrix(int8_t **mat, int64_t nrow, int64_t ncol) {
    int64_t i, j;
    for(i=0; i<nrow; i++) {
        for(j=0; j<ncol; j++) fprintf(stderr, "%"PRId8, mat[i][j]);
        fprintf(stderr, " %"PRId64"\n", i); fflush(stderr);
    }
}

void pushSGCIGAR(uint32_t **cigar, int32_t l, int32_t *m, int32_t op, int32_t oplen) {
    if(l+1 >= *m) {
        *m = l+1;
        kroundup32((*m));
        (*cigar) = realloc(*cigar, sizeof(uint32_t)*(*m));
        assert(*cigar);
    }
    (*cigar)[l] = bam_cigar_gen(oplen, op);
}

uint16_t getScore(int8_t *ref, int8_t *seq, int32_t seqLen, uint16_t maxScore, uint16_t scoreMatrix[4][4]) {
    int32_t i;
    uint16_t score = 0;
    for(i=0; i<seqLen; i++) {
        score += scoreMatrix[seq[i]][ref[i]];
        if(score >= maxScore) break;
    }
    return score;
}

//This method should average O(N*(M-N)/2 time, which is ~2x faster than GlobalAlignment().
//The likelyStartPos heuristic speeds this up significantly
s_align * GlobalAlignment(int8_t *ref, int32_t refLen, int8_t *path, int32_t pathLen, int k, int32_t likelyStartPos) {
    int32_t i, best_i = -1;
    uint16_t nmatch = 1, mismatch = 3, score;
    uint16_t best = mismatch * (pathLen>>2); //A heuristic to speed things up a bit
    uint32_t *cigar = malloc(sizeof(uint32_t));
    uint16_t scoreMatrix[4][4] = {{0, mismatch, mismatch, 0},
                                  {mismatch, 0, mismatch, 0},
                                  {mismatch, mismatch, 0, 0},
                                  {nmatch, nmatch, nmatch,0}};
    s_align *sal = malloc(sizeof(s_align));
    assert(sal);
    assert(cigar);
    cigar[0] = bam_cigar_gen(pathLen, 0);

    //Try the most likely start position
    score = getScore(ref+likelyStartPos, path, pathLen, best, scoreMatrix);
    if(score < best) {
        best = score;
        best_i = likelyStartPos;
    }

    //Try restricting the score range
    if(score != 0) {
        for(i=k; i<refLen-pathLen-k+1; i++) {
            if(i==likelyStartPos) continue;
            score = getScore(ref+i, path, pathLen, best, scoreMatrix);
            if(score < best) {
                best = score;
                best_i = i;
            }
        }
    }
    //If we didn't get a match, then use an unrestricted score range
    if(best_i == -1) {
        best = 0xFFFF;
        for(i=k; i<refLen-pathLen-k+1; i++) {
            score = getScore(ref+i, path, pathLen, best, scoreMatrix);
            if(score < best) {
                best = score;
                best_i = i;
            }
        }
    }

    sal->score1 = best;
    sal->ref_begin1 = best_i;
    sal->ref_end1 = best_i+pathLen-1;
    sal->read_begin1 = 0;
    sal->read_end1 = pathLen-1;
    sal->cigarLen = 1;
    sal->cigar = cigar;

    return sal;
}

//Perform semi-global alignment using a Needleman-Wunsch like algorithm. Since the first and last k of the 
//path and reference are known to be identical, this can ensure complete overlap. In the future, the
//mismatch, gap open, and gap extend penalties could be user specified.
//This could be made slightly faster by removing the topmost and leftmost k rows/columns, making the matrix (n-2k+1)*(m-2k+1)
s_align * SemiGlobalAlignment(int8_t *ref, int32_t refLen, int8_t *path, int32_t pathLen, int32_t k) {
    int64_t i, j, *last, *cur, *tmp;
    int8_t **mat;
    int64_t nmatch = -1, mismatch = -3, gapOpen = -5, gapExtend = -3, left, top, diag, best;
    int32_t op = 0, oplen = k, maxCigar = 1;
    uint32_t *cigar = NULL, tmpCIGAR;
    int64_t scoreMatrix[4][4] = {{0, mismatch, mismatch, 0},
                                 {mismatch, 0, mismatch, 0},
                                 {mismatch, mismatch, 0, 0},
                                 {nmatch, nmatch, nmatch,0}};
    s_align *sal = malloc(sizeof(s_align));
    assert(sal);

    sal->score1 = 0;
    sal->ref_begin1 = 0;
    sal->ref_end1 = refLen-1;
    sal->read_begin1 = 0;
    sal->read_end1 = pathLen-1;
    sal->cigarLen = 0;

    //If the sequences are identical, then life is easy
    if(refLen == pathLen) {
        for(i=0; i<pathLen;i++) if(ref[i] != path[i]) goto skip;
        sal->cigar = malloc(sizeof(uint32_t) * maxCigar);
        sal->cigar[0] = bam_cigar_gen(refLen, 0);
        sal->cigarLen = 1;
        return sal;
    }

skip:
    //Initialize the matrix and counts
    mat = malloc((pathLen+2-k)*sizeof(int8_t*));
    assert(mat);
    for(i=0; i<pathLen+2-k; i++) {
        mat[i] = calloc(refLen+2-k, sizeof(int8_t));
        assert(mat[i]);
    }
    for(j=1; j<refLen+2-k; j++) mat[0][j] = 4;
    last = malloc(sizeof(int64_t) * (refLen+2-k));
    cur = malloc(sizeof(int64_t) * (refLen+2-k));
    assert(last); assert(cur);
    last[0] = 0;
    for(j=1; j<refLen+2-k; j++) last[j] = gapExtend*j + gapOpen;

    //Fill in the direction matrix
    for(i=1; i<pathLen+2-k; i++) {
        mat[i][0] = 1;
        cur[0] = gapExtend*i + gapOpen;
        for(j=1; j<refLen+2-k; j++) {
            left = cur[j-1] + gapExtend + ((mat[i][j-1] & 4)?0:gapOpen);
            top = last[j] + gapExtend + ((mat[i-1][j] & 1)?0:gapOpen);
            diag = last[j-1] + scoreMatrix[path[i-1]][ref[j-1]];
            best = left;
            if(best < top) best = top;
            if(best < diag) best = diag;
            cur[j] = best;
            if(best == top) mat[i][j] |= 1;
            if(best == diag) mat[i][j] |= 2;
            if(best == left) mat[i][j] |= 4;
        }
        tmp = cur;
        cur = last;
        last = tmp;
    }

    //Trace-back
    cigar = malloc(sizeof(uint32_t) * maxCigar);
    assert(cigar);
    i = pathLen-k+1;
    j = refLen-k+1;
    oplen=k-1;
    while(i>k && j>k) {
        //Try to perform the last CIGAR op again
        if(op == 0) { //M, diag
            if(mat[i][j] & 2) {
                oplen++;
                i--;
                j--;
                continue;
            }
        } else if(op == 1) { //I, top
            if(mat[i][j] & 1) {
                oplen++;
                i--;
                continue;
            }
        } else { //D, left
            if(mat[i][j] & 4) {
                oplen++;
                j--;
                continue;
            }
        }
        //We've changed the op
        pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, op, oplen);
        oplen = 1;
        if(mat[i][j]&2) { //M
            op = 0;
            i--; j--;
        } else if(mat[i][j]&1) { //I
            op = 1;
            i--;
        } else { //D
            op = 2;
            j--;
        }
    }

    //Finish up
    if(i==k && j==k) {
        if(op != 0) {
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, op, oplen);
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, 0, k);
        } else {
            oplen += k;
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, op, oplen);
        }
    } else if(i==k) { //Still have a deletion
        if(op != 2) {
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, op, oplen);
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, 2, j-k);
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, 0, k);
        } else {
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, op, oplen+j-k);
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, 0, k);
        }
    } else { //Still have an insertion
        if(op != 1) {
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, op, oplen);
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, 1, i-k);
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, 0, k);
        } else {
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, op, oplen+i-k);
            pushSGCIGAR(&cigar, sal->cigarLen++, &maxCigar, 0, k);
        }
    }

    //Reverse the cigar order
    //It might be faster to malloc() and then iterate a fill, we should time both...
    for(i=0, j=sal->cigarLen-1; i<(sal->cigarLen)>>1; i++) {
        tmpCIGAR = cigar[i];
        cigar[i] = cigar[j];
        cigar[j] = tmpCIGAR;
        j--;
    }
    sal->cigar = cigar;
#ifdef DEBUG
//    printMatrix(mat, pathLen+2-k, refLen+2-k);
    fprintf(stderr, "[SemiGlobalAlignment] ref ");
    for(i=0; i<refLen; i++) fprintf(stderr, "%"PRId8, ref[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "[SemiGlobalAlignment] path ");
    for(i=0; i<pathLen; i++) fprintf(stderr, "%"PRId8, path[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "[SemiGlobalAlignment] CIGAR ");
    for(i=0; i<sal->cigarLen; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(sal->cigar[i]), BAM_CIGAR_STR[bam_cigar_op(sal->cigar[i])]);
    fprintf(stderr, "\n"); fflush(stderr);
    int32_t l = 0;
    for(i=0; i<sal->cigarLen; i++) {
        op = bam_cigar_op(sal->cigar[i]);
        oplen = bam_cigar_oplen(sal->cigar[i]);
        if(bam_cigar_type(op)&1) l+=oplen;
    }
    assert(l==pathLen);
    l = 0;
    for(i=0; i<sal->cigarLen; i++) {
        op = bam_cigar_op(sal->cigar[i]);
        oplen = bam_cigar_oplen(sal->cigar[i]);
        if(bam_cigar_type(op)&2) l+=oplen;
    }
    assert(l==refLen);
#endif

    for(i=0; i<pathLen+2-k; i++) free(mat[i]);
    free(mat);
    free(cur);
    free(last);
    return sal;
}
