#include <zlib.h>
#include "realigner.h"
#include "htslib/kseq.h"
#include "htslib/faidx.h"
#include "SSW/ssw.h"

KSTREAM_INIT(gzFile, gzread, 4096)

char int2base[32] = {0, 'A', 'T', 0, 'G', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N', \
                     0, 'A', 'C', 0, 'A', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N'};

                        // A CG  T  N
int8_t scoreMatrix[16] = { 2,-2,-2, 0,
                          -2, 2,-2, 0,
                          -2,-2, 2, 0,
                           0, 0, 0, 0};

struct bounds {
    int32_t tid, start, end;
};

//This is basically calculate_positions() from bison
//The output must be free()d
int32_t * bam2positions(bam1_t *b) {
    int32_t *positions, i = 0, j = 0;
    positions = malloc(sizeof(int32_t) * b->core.l_qseq);
    uint32_t *cigar = bam_get_cigar(b);
    int32_t prev = b->core.pos;
    int op, op_len, offset=0;

    for(i=0; i<b->core.n_cigar; i++) {
        op=bam_cigar_op(cigar[i]);
        op_len = bam_cigar_oplen(cigar[i]);
        switch(op) {
        case 0 : //M
        case 7 : //=
        case 8 : //X
            for(j=0; j<op_len; j++) {
                *(positions+offset) = prev++;
                offset++;
            }
            break;
        case 1 : //I
        case 4 : //S
            //Assign the previous position
            for(j=0; j<op_len; j++) {
                *(positions+offset) = prev;
                offset++;
            }
            break;
        case 2 : //D
        case 3 : //N
            prev += op_len;
            break;
        case 5 : //H
        case 6 : //P
            //Ignore these
            break;
        case 9 : //B
            //We don't validate this!
            prev -= op_len;
            break;
        default :
            fprintf(stderr, "[bam2positions] Ignoring an unknown CIGAR operation: %i\n", op);
            fflush(stderr);
            break;
        }
    }
    return positions;
}

//find the bounds of a read within k on either side of an ROI
void findPositions(int32_t *positions, int len, int32_t start, int32_t end, int *start2, int *end2) {
    int i = 0;
    *start2 = -1;
    for(i=0; i<len; i++) {
        if(positions[i] >= start && positions[i] <= end) {
            if(*start2 == -1) *start2 = i;
            *end2 = i;
        }
    }
}

//This is from PileOMeth
int getStrand(bam1_t *b) {
    char *XG = (char *) bam_aux_get(b, "XG");
    if(XG == NULL) { //Can't handle non-directional libraries!
        if(b->core.flag & BAM_FPAIRED) {
            if((b->core.flag & 0x50) == 0x50) return 2; //Read1, reverse comp. == OB
            else if(b->core.flag & 0x40) return 1; //Read1, forward == OT
            else if((b->core.flag & 0x90) == 0x90) return 1; //Read2, reverse comp. == OT
            else if(b->core.flag & 0x80) return 2; //Read2, forward == OB
            return 0; //One of the above should be set!
        } else {
            if(b->core.flag & 0x10) return 2; //Reverse comp. == OB
            return 0;
        }
    } else {
        if(*(XG+1) == 'C') { //OT or CTOT, due to C->T converted genome
            if((b->core.flag & 0x51) == 0x41) return 1; //Read#1 forward == OT
            else if((b->core.flag & 0x51) == 0x51) return 3; //Read #1 reverse == CTOT
            else if((b->core.flag & 0x91) == 0x81) return 3; //Read #2 forward == CTOT
            else if((b->core.flag & 0x91) == 0x91) return 1; //Read #2 reverse == OT
            else if(b->core.flag & 0x10) return 3; //Single-end reverse == CTOT
            else return 1; //Single-end forward == OT
        } else {
            if((b->core.flag & 0x51) == 0x41) return 4; //Read#1 forward == CTOB
            else if((b->core.flag & 0x51) == 0x51) return 2; //Read #1 reverse == OB
            else if((b->core.flag & 0x91) == 0x81) return 2; //Read #2 forward == OB
            else if((b->core.flag & 0x91) == 0x91) return 4; //Read #2 reverse == CTOB
            else if(b->core.flag & 0x10) return 2; //Single-end reverse == OB
            else return 4; //Single-end forward == CTOB
        }
    }
}

void bam2kmer(bam1_t *b, int k, int32_t start, int32_t end, bf *bf, kstring_t *ks, char *CT, char *GA, int32_t refLen) {
    int i, start2, end2;
    int offset = (getStrand(b) & 1) ? 0 : 16;
    int32_t *positions; //Could reuse this like a kstring_t
    int32_t refStart = (start-k>=0) ? start-k : 0;
    int32_t refEnd = refStart+refLen-1;
    uint64_t h;
#ifdef DEBUG
    char *tmp = calloc(k+1, sizeof(char));
#endif

    positions = bam2positions(b);
    findPositions(positions, b->core.l_qseq, start, end, &start2, &end2);

    //Grow ks as needed
    if(ks->m < (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart)+1)
        ks_resize(ks, (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart)+1);
    ks->l = (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart)+1;

    //Add 5' reference sequence
    if(positions[start2]-refStart) {
        if(offset==0) strncpy(ks->s, CT, positions[start2]-refStart);
        else strncpy(ks->s, GA, positions[start2]-refStart);
    }

    //Extract the C->T and G->A sequences
    for(i=start2; i<=end2; i++)
        ks->s[positions[start2]-refStart+i-start2] = int2base[offset+bam_seqi(bam_get_seq(b), i)];

    //Add 3' reference sequence
    if(refEnd-positions[end2]) {
        if(offset==0) strncpy(ks->s+positions[start2]-refStart+end2-start2+1,
                              CT+refLen-refEnd+positions[end2],
                              refEnd-positions[end2]);
        else strncpy(ks->s+positions[start2]-refStart+end2-start2+1,
                     GA+refLen-refEnd+positions[end2],
                     refEnd-positions[end2]);
    }
    ks->s[ks->l] = '\0';
#ifdef DEBUG
    fprintf(stderr, "[bam2kmer] Final read sequence %s\n", ks->s); fflush(stderr);
#endif

    //Add the kmers
    for(i=0; i<ks->l-k; i++) {
        h = hash_seq(ks->s+i, k);
#ifdef DEBUG
        snprintf(tmp, k+1, "%s", ks->s+i);
        fprintf(stderr, "[bam2kmer] vertex %s %" PRIu64"\n", tmp, h & bf->mask);
#endif
        bf_add(bf, h);
    }

#ifdef DEBUG
    free(tmp);
#endif
    free(positions);
}

void convertCT(char *s, int l) {
    int i;
    for(i=0; i<l; i++) {
        switch(s[i]) {
        case 'A' :
        case 'a' :
            s[i] = 'A';
            break;
        case 'C' :
        case 'c' :
        case 'T' :
        case 't' :
            s[i] = 'T';
            break;
        case 'G' :
        case 'g' :
            s[i] = 'G';
            break;
        default :
            s[i] = 'N';
            break;
        }
    }
}

void convertGA(char *s, int l) {
    int i;
    for(i=0; i<l; i++) {
        switch(s[i]) {
        case 'A' :
        case 'a' :
        case 'G' :
        case 'g' :
            s[i] = 'A';
            break;
        case 'C' :
        case 'c' :
            s[i] = 'C';
            break;
        case 'T' :
        case 't' :
            s[i] = 'T';
            break;
        default :
            s[i] = 'N';
            break;
        }
    }
}

int8_t *char2int(char *seq, int len) {
    int i;
    int8_t *o = malloc(sizeof(int8_t) * len);
    assert(o);

    for(i=0; i<len; i++) {
        switch(seq[i]) {
        case 'A' :
            o[i] = 0;
            break;
        case 'C' :
        case 'G' :
            o[i] = 1;
            break;
        case 'T' :
            o[i] = 2;
            break;
        case 'N' :
            o[i] = 3;
            break;
        default :
            fprintf(stderr, "[char2int] Replacing unknown base '%c' with 'N'!\n", seq[i]);
            fflush(stderr);
            o[i] = 3;
            break;
        }
    }
    return o;
}

//Given a set of alignments and a score, return the alignment with that score with the highest count
int32_t findBestAlignment(s_align **al, int32_t len, uint32_t *counts) {
    int32_t best=-1, bestScore = -1, i;
    uint32_t bestCount = 0;

    //Get the maximum score
    for(i=0; i<len; i++)
        if(al[i]->score1 > bestScore && al[i]->cigarLen == 1) bestScore = al[i]->score1;

    //Given a score, find the path with the highest count
    for(i=0; i<len; i++) {
        if(al[i]->score1 == bestScore && counts[i] > bestCount) {
            bestCount = counts[i];
            best = i;
        }
    }

    return best;
}

//Count the number of best-stratum alignments per path
void countAlignmentsPerPath(s_align** al, int32_t len, uint32_t *counts) {
    int32_t bestScore = -1, i;

    //Get the maximum score
    for(i=0; i<len; i++)
        if(al[i]->score1 > bestScore && al[i]->cigarLen == 1) bestScore = al[i]->score1;

#ifdef DEBUG
    fprintf(stderr, "[countAlignmentsPerPath] bestScore is %" PRId32 "\n", bestScore);
#endif
    //Increment the counter if score1==bestScore
    for(i=0; i<len; i++) {
        if(al[i]->score1 == bestScore && al[i]->cigarLen == 1) {
            if(counts[i] < 0xFFFFFFFF) counts[i]++;
#ifdef DEBUG
            fprintf(stderr, "[countAlignmentsPerPath] Incrementing %" PRId32 "\n", i);
#endif
        }
    }
}

//Update b->core.pos, adding an OP aux tag if needed
void updatePos(bam1_t *b, int32_t newStartPos) {
    uint8_t *p;
    int32_t oldpos = 1+b->core.pos; //1-based
    if(b->core.pos != newStartPos) {
        p = bam_aux_get(b, "OP");
        if(p==NULL) bam_aux_append(b, "OP", 'i', 4, (uint8_t*) &oldpos);
        b->core.pos = newStartPos;
    }
}

//Update the CIGAR string, appending an OC if needed
void updateCigar(bam1_t *b, int32_t *arr, int len) {
    uint8_t *p;
    size_t remaining_data_len;
    int32_t i;
    int addOC=0;

    if(len == b->core.n_cigar) {
        for(i=0; i<len; i++) {
            if(arr[i] != bam_get_cigar(b)[i]) {
                addOC=1;
                break;
            }
        }
    } else addOC=1;

    if(addOC) {
        p = bam_aux_get(b, "OC");
        if(p==NULL) {
            uint32_t *cigar = bam_get_cigar(b);
            kstring_t *str = calloc(1, sizeof(kstring_t));
            for (i = 0; i < b->core.n_cigar; ++i) {
                kputw(bam_cigar_oplen(cigar[i]), str);
                kputc(bam_cigar_opchr(cigar[i]), str);
            }
            bam_aux_append(b, "OC", 'Z', str->l+1, (uint8_t*) str->s);
        }
    }

    //Do we need to move everything after the cigar?
    if(b->core.n_cigar < len) {
        //Do we need to realloc b->data?
        if(b->m_data-b->l_data > 4*(len-b->core.n_cigar)) {
            b->m_data++;
            kroundup32(b->m_data);
            b->data = realloc(b->data, b->m_data*sizeof(uint8_t));
            assert(b->data);
        }
        p = bam_get_seq(b);
        remaining_data_len = b->data+b->l_data-p;
        memmove(p+4*(len-b->core.n_cigar), p, remaining_data_len);
        b->l_data += 4*(len-b->core.n_cigar);
    }
    for(i=0; i<len; i++) ((int32_t *) bam_get_cigar(b))[i] = arr[i];
    b->core.n_cigar = len;
}

int32_t *pushCIGAR(int32_t *CIGAR, int32_t n_cigar, int32_t *max_cigar, int32_t op, int32_t oplen) {
    if(n_cigar+1 >= *max_cigar) {
        (*max_cigar)++;
        kroundup32((*max_cigar));
        CIGAR= realloc(CIGAR, (*max_cigar) * sizeof(int32_t));
        assert(CIGAR);
    }
    CIGAR[n_cigar] = bam_cigar_gen(oplen, op);
    return CIGAR;
}
    
/*
readStartPos	0-based position in the read that denotes the 5'-most base used in the alignment
readEndPos	0-based position in the read that denotes the 3'-most base used in the alignment
refStart	The 0-based genomic coordinate of the 5'-most base of the reference pathway
refStartPos	0-based position in the reference that denotes the 5'-most base used in the alignment
refEndPos	0-based position in the reference that denotes the 3'-most base used in the alignment
*/
bam1_t * updateAlignment(bam1_t *b, s_align *al, int32_t readStartPos, int32_t readEndPos, int32_t refStartPos, int32_t refStart, int32_t refEndPos) {
    int32_t *newCIGARArray = malloc(sizeof(uint32_t) * b->core.n_cigar);
    int n_cigar = 0, max_cigar = b->core.n_cigar, read_cigar_opnum = 0, ref_cigar_opnum = 0;
    int32_t oplen = 0, op,oplen2=0, op2=0, readPos = 0, newStartPos = -1;
    int32_t i, refPos = 0;

    assert(newCIGARArray);

    /* We deal with things in 4 segments:
     * (1) read 5' overhang
     * (2) read 5' underhang
     * (3) overlap region
     * (4) read 3' overhang
     * Note that a read 3' underhang can be ignored!
     */
//    fprintf(stderr, "readStartPos %" PRId32 " readpos %" PRId32 " readEndPos %" PRId32 "\n", readStartPos, readPos, readEndPos);
//    fprintf(stderr, "refStartPos %" PRId32 "\n", refStartPos);
    if(readStartPos > 0) { //5' overhang
        while(readPos < readStartPos) {
            if(oplen == 0) {
                op = bam_cigar_op(bam_get_cigar(b)[read_cigar_opnum]);
                oplen = bam_cigar_oplen(bam_get_cigar(b)[read_cigar_opnum++]);
            }
            //If the op consumes the query and extends past the while-loop bounds, then trim
//            fprintf(stderr, "readPos %" PRId32 " oplen %"PRId32" readStartPos %"PRId32" bam_cigar_type(op) %i\n", readPos, oplen, readStartPos, bam_cigar_type(op));
            if(readPos + oplen>= readStartPos && (bam_cigar_type(op)&1)) {
                oplen = readStartPos-readPos;
                break;
            } else { //Append to newCIGARArray
                newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
                oplen = 0;
            }
        }
        //If refStartPos >0, we need to cat up to it
        if(refStartPos > 0) {
            while(refPos < refStartPos) {
                if(oplen2 == 0) {
                    op2 = bam_cigar_op(al->cigar[ref_cigar_opnum]);
                    oplen2 = bam_cigar_oplen(al->cigar[ref_cigar_opnum++]);
                }
                if(refPos + oplen2>= refStartPos && bam_cigar_type(op)&2) {
                    oplen2 = refPos+oplen2-refStartPos;
                    refPos = refStartPos;
                    break;
                }
            }
        }
//    fprintf(stderr, "A CIGAR ");
//    for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
//    fprintf(stderr, "\n");
    } else if(refStartPos>0) { //5' Underhang, note new start position
        newStartPos = refStart;
//        fprintf(stderr, "refStart is %" PRId32 " so we start with newStartPos = that\n", refStart);
        while(refPos < refStartPos) {
            if(oplen2 == 0) {
                op2 = bam_cigar_op(al->cigar[ref_cigar_opnum]);
                oplen2 = bam_cigar_oplen(al->cigar[ref_cigar_opnum++]);
            }
//            fprintf(stderr, "refPos %" PRId32 " refStartPos %" PRId32 "\n", refPos, refStartPos);
//            fprintf(stderr, "oplen %" PRId32 " op %c\n", oplen2, BAM_CIGAR_STR[op2]);
//            fprintf(stderr, "bam_cigar_type(op) == %i\n", bam_cigar_type(op2));
//            fflush(stderr);
            if(refPos + oplen2>= refStartPos && bam_cigar_type(op2)&1) {
                if(bam_cigar_type(op2) &2) newStartPos += refStartPos;
                oplen2 = refPos+oplen2-refStartPos;
//                fprintf(stderr, "Truncating oplen to %" PRId32 "\n", oplen2);
                refPos = refStartPos;
                break;
            }
            if(bam_cigar_type(op2) &2) {//OP consumes the reference
                newStartPos += oplen2;
            }
            if(bam_cigar_type(op2) &1) {//OP consumes the query
                refPos += oplen2;
            }
            oplen2 = 0;
        }
//    fprintf(stderr, "B CIGAR ");
//    for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
//    fprintf(stderr, "\n");
    }
//    fprintf(stderr, "Remnant CIGAR is %" PRId32 "%c\n", oplen, BAM_CIGAR_STR[op]);
//    fprintf(stderr, "Remnant CIGAR2 is %" PRId32 "%c\n", oplen2, BAM_CIGAR_STR[op2]);
//    fflush(stderr);

    //Deal with the overlap region
//    fprintf(stderr, "refPos %" PRId32 " refEndPos %" PRId32 "\n", refPos, refEndPos);
    while(refPos <= refEndPos) {
        if(oplen2 == 0) {
            op2 = bam_cigar_op(al->cigar[ref_cigar_opnum]);
            oplen2 = bam_cigar_oplen(al->cigar[ref_cigar_opnum++]);
        }
//        fprintf(stderr, "1 Found op %"PRId32"%c bam_cigar_type(op2) %i\n", oplen2, BAM_CIGAR_STR[op2], bam_cigar_type(op2));
//        fprintf(stderr, "1 refPos %" PRId32 " refEndPos+1 %" PRId32 "\n", refPos, refEndPos+1);
        if(refPos + oplen2>= refEndPos+1 && bam_cigar_type(op2)&1) {
            oplen2 = refEndPos-refPos+1;
//            fprintf(stderr, "Truncating the operation to length %" PRId32" and breaking\n", oplen2);
            if(oplen) {
                if(op == op2) {
                    oplen2 += oplen;
                } else newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
                oplen = 0;
            }
            refPos += oplen2;
            readPos += oplen2;
            break;
        } else { //Append to newCIGARArray
            if(oplen) {
//                fprintf(stderr, "Comparing %c and %c\n", BAM_CIGAR_STR[op], BAM_CIGAR_STR[op2]);
                if(op == op2) {
//                    fprintf(stderr, "refPos before %" PRId32" ", refPos);
                    if(bam_cigar_type(op2) &1) {
                        refPos += oplen2;
                        readPos += oplen2;
                    }
//                    fprintf(stderr, " and %" PRId32" after\n", refPos);
                    oplen2 += oplen;
                } else newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
                oplen = 0;
            } else {
                if(bam_cigar_type(op2) &1) refPos += oplen2;
            }
//            fprintf(stderr, "pushing %"PRId32"%c\n", oplen2, BAM_CIGAR_STR[op2]);
            newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op2, oplen2);
            oplen2 = 0;
        }
        oplen2 = 0;
//        fprintf(stderr, "2 refPos %" PRId32 " refEndPos+1 %" PRId32 "\n", refPos, refEndPos+1);
    }
//    fprintf(stderr, "Remnant CIGAR is %" PRId32 "%c\n", oplen2, BAM_CIGAR_STR[op2]);
    oplen = 0;
//    fprintf(stderr, "C CIGAR ");
//    for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
//    fprintf(stderr, "\n");

    //3' Overhang
//    fprintf(stderr, "Before D, readEndPos %" PRId32 " b->core.l_qseq-1 %" PRId32"\n", readEndPos, b->core.l_qseq-1);
    if(readEndPos < b->core.l_qseq-1) {
        //Iterate up to the relevant CIGAR operation
        i = 0; //within-read position index
        for(read_cigar_opnum=0; read_cigar_opnum<b->core.n_cigar; read_cigar_opnum++) {
            op = bam_cigar_op(bam_get_cigar(b)[read_cigar_opnum]);
            oplen = bam_cigar_oplen(bam_get_cigar(b)[read_cigar_opnum]);
            if(bam_cigar_type(op) & 1) {
                if(i+oplen>=readEndPos) {
                    oplen = oplen+i-readEndPos;
                    readPos += oplen;
                    break;
                }
                readPos += oplen;
            }
        }
//        fprintf(stderr, "D CIGAR ");
//        for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
//        fprintf(stderr, "\n");
        //Can we merge this OP with that from the overlap?
        if(oplen2) {
            if(op==op2) {
                oplen +=oplen2;
            } else {
                newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op2, oplen2);
                oplen2 = 0;
            }
        }
//        fprintf(stderr, "E CIGAR ");
//        for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
//        fprintf(stderr, "\n");
        //Finish adding the alignment's original OPs
        for(; read_cigar_opnum<b->core.n_cigar; read_cigar_opnum++) {
            op = bam_cigar_op(bam_get_cigar(b)[read_cigar_opnum]);
            oplen = bam_cigar_oplen(bam_get_cigar(b)[read_cigar_opnum]);
            newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
        }
    } else if(oplen2) { //Push the last overlap CIGAR op
        newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op2, oplen2);
    }

#ifdef DEBUG
    fprintf(stderr, "[updateAlignment] Final CIGAR ");
    for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
    fprintf(stderr, "\n");
    if(newStartPos != -1) fprintf(stderr, "[updateAlignment] new alignment starting position %" PRId32 "\n", newStartPos);
    fflush(stderr);
#endif
    if(newStartPos != -1) updatePos(b, newStartPos);
    updateCigar(b, newCIGARArray, n_cigar);
    return b;
}

//Align all of the paths to a reference sequence
s_align **alignPaths2Ref(paths *p, int8_t *ref, int32_t refLen, char *refSeq, int32_t *refIndex) {
    int32_t i;
    s_profile *s;
    s_align **alignments = malloc(sizeof(s_align*) * p->l);
    assert(alignments);
#ifdef DEBUG
    int j;
#endif
    for(i=0; i<p->l; i++) {
        s = ssw_init(p->conv[i], p->len[i], scoreMatrix, 4, (2*p->len[i] < 255) ? 0: 1);
        alignments[i] = ssw_align(s, ref, refLen, 3, 1, 2, 0, 0, 15);
#ifdef DEBUG
        fprintf(stderr, "[alignPaths2Ref] path[%i] aligned against ref with CIGAR ", i);
        for(j=0; j<alignments[i]->cigarLen; j++) fprintf(stderr, "%" PRIu32 "%c", cigar_int_to_len(alignments[i]->cigar[j]), cigar_int_to_op(alignments[i]->cigar[j]));
        fprintf(stderr, "\n");
        fflush(stderr);
#endif
        if(*refIndex == -1 && alignments[i]->cigarLen == 1)
            if(strcmp(refSeq, p->path[i]) == 0) *refIndex = i;
        init_destroy(s);
    }
    return alignments;
}

//Align a read to the appropriate set of paths
s_align **alignReads2Paths(bam1_t *b, int strand, int32_t *subreadM, int8_t **subreadSeq, paths *p, int32_t refLBound, int32_t refRBound, int32_t *readLBound, int32_t *readRBound) {
                           //   A  C     G           T                    N
    int8_t int2small[16] = {-1, 0, 1,-1, 1,-1,-1, 0, 2,-1,-1,-1,-1,-1,-1, 3};
    int32_t *positions, i, subreadL;
    s_profile *s;
    s_align **readal = malloc(sizeof(s_align*)*p->l);
    assert(readal);
#ifdef DEBUG
    int j;
#endif

    //Make int2small AGTN or ACTN
    int2small[2] = (strand&1) ? 2 : 1; //C->A
    int2small[4] = ((strand&1) == 0) ? 0 : 1; //G->A

    //extract int8_t * covering the ref
    positions = bam2positions(b);
    findPositions(positions, b->core.l_qseq, refLBound, refRBound, readLBound, readRBound);
    subreadL = (*readRBound)-(*readLBound)+1;
    if(subreadL >= *subreadM) { //Grow subreadSeq if needed
        (*subreadM)++;
        kroundup32((*subreadM));
        (*subreadSeq) = realloc((*subreadSeq), sizeof(int8_t) * (*subreadM));
        assert((*subreadSeq));
    }
    for(i=(*readLBound); i<=(*readRBound); i++) (*subreadSeq)[i-(*readLBound)] = int2small[bam_seqi(bam_get_seq(b), i)];
#ifdef DEBUG

    fprintf(stderr, "[alignReads2Paths] read ");
    for(i=0; i<subreadL; i++) fprintf(stderr, "%"PRId8, (*subreadSeq)[i]);
    fprintf(stderr, "\n");
    fflush(stderr);
#endif

    //Align
    s = ssw_init((*subreadSeq), subreadL, scoreMatrix, 4, (2*subreadL < 255) ? 0: 1);
    for(i=0; i<p->l; i++) {
        readal[i] = ssw_align(s, p->conv[i], p->len[i], 3, 1, 2, 0, 0, 15);
#ifdef DEBUG
        fprintf(stderr, "[alignReads2Paths] ref[%i] ", i);
        for(j=0; j<p->len[i]; j++) fprintf(stderr, "%"PRId8, p->conv[i][j]);
        fprintf(stderr, "\n");
        fprintf(stderr, "[alignReads2Paths] readal[%i]->score1 %" PRIu16 "\n", i, readal[i]->score1);
        fprintf(stderr, "[alignReads2Paths] readal[%i]->ref_begin1 %" PRId32 " ref_end1 %" PRId32, i, readal[i]->ref_begin1, readal[i]->ref_end1);
        fprintf(stderr, "[alignReads2Paths] read_begin1 %" PRId32 " read_end1 %" PRId32 "\n", readal[i]->read_begin1, readal[i]->read_end1);
        fprintf(stderr, "[alignReads2Paths] CIGAR ");
        for(j=0; j<readal[i]->cigarLen; j++) fprintf(stderr, "%" PRIu32 "%c", cigar_int_to_len(readal[i]->cigar[j]), cigar_int_to_op(readal[i]->cigar[j]));
        fprintf(stderr, "\n");
        fflush(stderr);
#endif
    }

    //Clean up
    free(positions);
    init_destroy(s);

    return(readal);
}

//This should be made multithreaded);
void realignHeapCore(alignmentHeap **heap, paths *CTpaths, paths *GApaths, char *CT, char *GA, int k) {
    int32_t i, refLen = strlen(CT), refIndex[2] = {-1, -1};
    int8_t *CTref, *GAref, *subreadSeq;
    s_align **CTal, **GAal, ***readal;
    int32_t refLBound, refRBound;
    int32_t subreadM, bestAl;
    int strand;
    int32_t *readLBound = malloc(sizeof(int32_t) * (*heap)->l);
    int32_t *readRBound = malloc(sizeof(int32_t) * (*heap)->l);
    assert(readLBound);
    assert(readRBound);

    //subreadSeq and subreadM combine to function similar to a kstring_t
    //we don't need to constantly malloc() space to hold part of a read (subreadSeq)
    subreadSeq = malloc(50*sizeof(int8_t));
    assert(subreadSeq);
    subreadM = 50;

    readal = malloc((*heap)->l * sizeof(s_align**));
    assert(readal);

    //setup the sequences
    CTref = char2int(CT, refLen);
    GAref = char2int(GA, refLen);
    CTpaths->conv = malloc(sizeof(int8_t*)*CTpaths->l);
    GApaths->conv = malloc(sizeof(int8_t*)*GApaths->l);
    assert(CTpaths->conv);
    assert(GApaths->conv);
    for(i=0; i<CTpaths->l; i++) CTpaths->conv[i] = char2int(CTpaths->path[i], CTpaths->len[i]);
    for(i=0; i<GApaths->l; i++) GApaths->conv[i] = char2int(GApaths->path[i], GApaths->len[i]);

    //Align the paths to the reference sequence
    CTal = alignPaths2Ref(CTpaths, CTref, refLen, CT, &(refIndex[0]));
    GAal = alignPaths2Ref(GApaths, GAref, refLen, GA, &(refIndex[1]));

    //Iterate over the reads, aligning them to the paths
    refLBound = ((*heap)->start-k >= 0) ? (*heap)->start-k : 0;
    refRBound = refLBound + refLen-2*k;
    for(i=0; i<(*heap)->l; i++) {
        strand = getStrand((*heap)->heap[i]);
        readal[i] = alignReads2Paths((*heap)->heap[i], strand, &subreadM, &subreadSeq, (strand==1) ? CTpaths : GApaths, refLBound, refRBound, readLBound+i, readRBound+i);
    }

    //Count number of best/equally best alignments/path
    uint32_t *CTcounts = calloc(CTpaths->l, sizeof(uint32_t));
    uint32_t *GAcounts = calloc(GApaths->l, sizeof(uint32_t));
    for(i=0; i<(*heap)->l; i++) {
        strand = getStrand((*heap)->heap[i]);
        if(strand==1) countAlignmentsPerPath(readal[i], CTpaths->l, CTcounts);
        else countAlignmentsPerPath(readal[i], GApaths->l, GAcounts);
    }
#ifdef DEBUG
    for(i=0;i<CTpaths->l; i++) {
        fprintf(stderr, "CTcounts[%i] %" PRIu32 "\n", i, CTcounts[i]);
    }
    for(i=0;i<GApaths->l; i++) {
        fprintf(stderr, "GAcounts[%i] %" PRIu32 "\n", i, GAcounts[i]);
    }
#endif

    //Pick a best path for each gene given all of the others and update read
    refLBound = ((*heap)->start-2*k >= 0) ? (*heap)->start-2*k : 0; //The previous bounds were only for extracting subreads
    for(i=0; i<(*heap)->l; i++) {
        strand = getStrand((*heap)->heap[i]);
        if(strand&1) {
            bestAl = findBestAlignment(readal[i], CTpaths->l, CTcounts);
        } else {
            bestAl = findBestAlignment(readal[i], GApaths->l, GAcounts);
        }
        //Update read if it aligned to an alternate path!
#ifdef DEBUG
        fprintf(stderr, "[realignHeapCore] bestAl is #%" PRId32 "\n", bestAl);
        fprintf(stderr, "[realignHeapCore] %s readLBound %" PRId32" readRBound %" PRId32 "\n", bam_get_qname((*heap)->heap[i]), readLBound[i], readRBound[i]);
        fflush(stderr);
#endif
        if(bestAl == -1) continue;
        if(strand&1 && refIndex[0] != bestAl)
            (*heap)->heap[i] = updateAlignment((*heap)->heap[i],
                                               CTal[bestAl],
                                               readLBound[i],
                                               readRBound[i],
                                               readal[i][bestAl]->ref_begin1,
                                               refLBound,
                                               readal[i][bestAl]->ref_end1);
        else if(refIndex[1] != bestAl)
            (*heap)->heap[i] = updateAlignment((*heap)->heap[i],
                                               GAal[bestAl],
                                               readLBound[i],
                                               readRBound[i],
                                               readal[i][bestAl]->ref_begin1,
                                               refLBound,
                                               readal[i][bestAl]->ref_end1);
    }

    free(CTcounts);
    free(GAcounts);
    free(CTref);
    free(GAref);
    int j;
    for(i=0; i<(*heap)->l; i++) {
        strand = getStrand((*heap)->heap[i]);
        if(strand&1) {
            for(j=0; j<CTpaths->l; j++) align_destroy(readal[i][j]);
        } else {
            for(j=0; j<GApaths->l; j++) align_destroy(readal[i][j]);
        }
        free(readal[i]);
    }
    free(readal);
    for(i=0; i<CTpaths->l; i++) align_destroy(CTal[i]);
    for(i=0; i<GApaths->l; i++) align_destroy(GAal[i]);
}

//k is the kmer
void realignHeap(alignmentHeap *heap, int k, faidx_t *fai) {
    bf *filter = bf_init(heap->end-heap->start, WINDOW);
    int32_t i, start, end, start2;
    kstring_t *ks = calloc(1, sizeof(kstring_t));
    char *CT, *GA, *startVertex;
    paths *CTpaths, *GApaths;
    int len;
    vertex *CTcycles, *GAcycles;
    uint64_t h;

    //Add the reference sequence
    start = (heap->start-k >= 0) ? heap->start-k : 0;
    end = heap->end+k;
    start2 = start-k;
    if(start2<0) start2 = 0;
    CT = faidx_fetch_seq(fai, faidx_iseq(fai, heap->heap[0]->core.tid), start2, end+k-1, &len);
    GA = strdup(CT);
    convertCT(CT, len);
    convertGA(GA, len);
#ifdef DEBUG
    char *tmp = calloc(k+1, sizeof(char));
    fprintf(stderr, "[realignHeap] C->T reference %s\n", CT);
    fprintf(stderr, "[realignHeap] G->A reference %s\n", GA);
#endif
    for(i=0; i<len-k+1; i++) {
#ifdef DEBUG
        snprintf(tmp, k+1, "%s", CT+i);
        fprintf(stderr, "[realignHeap] CTvertex %s\n", tmp); fflush(stderr);
#endif
        h = hash_seq(CT+i, k);
        bf_add(filter, h);
        h = hash_seq(GA+i, k);
#ifdef DEBUG
        snprintf(tmp, k+1, "%s", GA+i);
        fprintf(stderr, "[realignHeap] GAvertex %s\n", tmp); fflush(stderr);
#endif
        bf_add(filter, h);
    }

    //Process the reads
    for(i=0; i<heap->l; i++) bam2kmer(heap->heap[i], k, start, end, filter, ks, CT, GA, len);

    //Extract the graph paths
    startVertex = strndup(CT, k);
    CTcycles = getCycles(filter, startVertex, CT+len-k, k, 'G');
#ifdef DEBUG
    fprintf(stderr, "[realignHeap] Going from %s -> %s\n", startVertex, CT+len-k); fflush(stdout);
#endif
    CTpaths = getPaths(filter, startVertex, CT+len-k, &CTcycles, 'G');
    snprintf(startVertex, (k+1)*sizeof(char), "%s", GA);
    GAcycles = getCycles(filter, startVertex, GA+len-k, k, 'C');
#ifdef DEBUG
    fprintf(stderr, "[realignHeap] Going from %s -> %s\n", startVertex, GA+len-k); fflush(stdout);
#endif
    GApaths = getPaths(filter, startVertex, GA+len-k, &GAcycles, 'C');
    free(startVertex);

    //Realign the read portions to the paths
    if(CTpaths->l + GApaths->l) {
        realignHeapCore(&heap, CTpaths, GApaths, CT, GA, k);
    } else {
        fprintf(stderr, "[realignHeap] Skipping %s:%" PRId32 "-%" PRId32 ", couldn't find any paths post-assembly!\n", faidx_iseq(fai, heap->heap[0]->core.tid), heap->start, heap->end);
    }
    //Clean up
    free(ks->s);
    free(ks);
    if(CTpaths) destroyPaths(CTpaths);
    if(GApaths) destroyPaths(GApaths);
    if(CTcycles) destroyDFSLL(CTcycles);
    if(GAcycles) destroyDFSLL(GAcycles);
    bf_destroy(filter);
    free(CT);
    free(GA);
}

//Move everything below here to a new file




//Does this really work correctly? The bounds are 0-base semi-open...
int regionNodeCmp(int32_t tid, int32_t start, int32_t end) {
    if(tid == lastTargetNode->tid) {
        if((start <= lastTargetNode->start && end >= lastTargetNode->start) || \
            (start >= lastTargetNode->start && start <= lastTargetNode->end)) return 0;
        else if(end < lastTargetNode->start) return end-lastTargetNode->start;
        assert(start > lastTargetNode->end);
        return start - lastTargetNode->end;
    } else return tid - lastTargetNode->tid;
}

//Does this really work correctly? The bounds are 0-base semi-open...
int overlapsRegionsCore(int32_t tid, int32_t start, int32_t end) {
    int direction;

    //The only reliable way to find the 5'-most region is from the left.
    while(regionNodeCmp(tid,start,end) <= 0) {
#ifdef DEBUG
        fprintf(stderr, "[overlapsRegionsCore] Direction %i\n", regionNodeCmp(tid,start,end));
        fprintf(stderr, "[overlapsRegionsCore] CurPosition  %i:%i-%i\n", (int) tid, (int)start, (int)end);
        fprintf(stderr, "[overlapsRegionsCore] ListPosition %i:%i-%i\n", (int) lastTargetNode->tid, (int)lastTargetNode->start, (int)lastTargetNode->end);
        fflush(stderr);
#endif
        lastTargetNode = lastTargetNode->prev; //Not sure if this really updates things
    }
    while(lastTargetNode->next != NULL) {
        lastTargetNode = lastTargetNode->next;
        direction = regionNodeCmp(tid,start,end);
        if(direction == 0) return 1;
        else if(direction > 0) return 0;
    }
    return 0;
}

//Does this really work correctly? The bounds are 0-base semi-open...
struct bounds * overlapsRegion(bam1_t *b) {
    int32_t tid = b->core.tid;
    int32_t start = b->core.pos;
    int32_t end = bam_endpos(b)-1;
    uint32_t *cigar = bam_get_cigar(b);
    int i, op, oplen;
    struct bounds *bounds = NULL;

    //Adjust bounds in the case of soft-clipping
    for(i=0; i<b->core.n_cigar; i++) {
        op = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        if(op == 4) {
            start -= oplen;
        } else break;
    }
    if(start < 0) start = 0;

    if(overlapsRegionsCore(tid, start, end)) {
        bounds = malloc(sizeof(struct bounds));
        assert(bounds != NULL);
        bounds->tid = lastTargetNode->tid;
        bounds->start = lastTargetNode->start;
        bounds->end = lastTargetNode->end;
    }
    return bounds;
}

void processReads(htsFile *fp, bam_hdr_t *hdr, htsFile *of, int k, faidx_t *fai) {
    bam1_t *b = bam_init1();
    struct bounds *bounds;
    alignmentHeap *heap = alignmentHeap_init(1000); //Need to make this an option
    int status; //1: EOF, 2: heap max, 3: Past ROI
    InDel reg;

    while(sam_read1(fp, hdr, b) > 0) {
#ifdef DEBUG
        fprintf(stderr, "[processReads] Found %s\n", bam_get_qname(b)); fflush(stderr);
#endif
        if(!(b->core.flag&4) && (bounds = overlapsRegion(b)) != NULL) {
#ifdef DEBUG
            fprintf(stderr, "[processReads] %s in region\n", bam_get_qname(b)); fflush(stderr);
#endif
            heap->heap[0] = b;
            heap->start = bounds->start;
            heap->end = bounds->end;
            heap->l = 1;
            b = bam_init1();
            assert(b != NULL);
            free(bounds);

            //We're in an ROI, add to the heap until:
            // (1) It's too big and we just skip the alignment
            // (2) We're past the 3' bound of the ROI
            status = 0;
inheap:     while(1) {
                if(sam_read1(fp, hdr, b)<0) {
                    status = 1;
                    break;
                }
#ifdef DEBUG
                fprintf(stderr, "[processReads] Found %s while in heap of length %i\n", bam_get_qname(b), heap->l); fflush(stderr);
                fprintf(stderr, "b->core.pos %" PRId32 " headp->end %" PRId32 "\n", b->core.pos, heap->end);
#endif
                if(b->core.pos < heap->end && b->core.tid == heap->heap[0]->core.tid) {
                    //If the alignment and the first in the heap have the same pos...
                    if(b->core.pos == heap->heap[0]->core.pos) {
                        reg.tid = b->core.tid;
                        reg.start = b->core.pos;
                        reg.end = bam_endpos(b)+1;
                        if(TargetNodeCmp(&reg, lastTargetNode) != 0) {
                            sam_write1(of, hdr, b);
                            continue;
                        }
                    }
                    if(heap->l+1 > heap->m) {
                        status = 2;
                        break;
                    } else {
                        heap->heap[heap->l++] = b;
                        b = bam_init1();
                        assert(b != NULL);
                    }
                } else {
                    status = 3;
                    break;
                }
            }

            //Deal with the 3 possible statuses
            assert(status != 0);
            if(status == 1) {
#ifdef DEBUG
                fprintf(stderr, "[processReads] Reached end of alignments\n"); fflush(stderr);
#endif
                if(heap->l) realignHeap(heap, k, fai);
#ifdef DEBUG
                fprintf(stderr, "[processReads] Last heap realigned\n"); fflush(stderr);
#endif
                writeHeap(of, hdr, heap);
                heap->l = 0;
                break;
            } else if(status == 2) {
                fprintf(stderr, "[processReads] Skipping %s:%"PRId32"-%"PRId32", too many reads\n", hdr->target_name[heap->heap[0]->core.tid], heap->start, heap->end);
            } else {
#ifdef DEBUG
                fprintf(stderr, "[processReads] %s is outside the ROI, realigning the heap\n", bam_get_qname(b)); fflush(stderr);
#endif
                realignHeap(heap, k, fai);
            }
            lastTargetNode = lastTargetNode->next; //Move to the next ROI
            heap = writeHeapUntil(of, hdr, heap);
            if(heap->l) {//The heap wasn't flushed
                goto inheap; //Yeah yeah, I know
            }
#ifdef DEBUG
        } else {
            fprintf(stderr, "[processReads] %s not in ROI\n", bam_get_qname(b)); fflush(stderr);
#endif
        }
        sam_write1(of, hdr, b);
    }
    bam_destroy1(b);
    alignmentHeap_destroy(heap);
}

//return -1 on error
int bed2list(gzFile fp, bam_hdr_t *hdr) {
    int ret;
    char *s;
    InDel *node;
    kstream_t *kstr = ks_init(fp);
    kstring_t *ks = calloc(1, sizeof(kstring_t));

    if(ks == NULL) return -1;
    lastTargetNode = firstTargetNode;
    while(ks_getuntil(kstr, '\n', ks, &ret) >= 0) {
        node = calloc(1, sizeof(InDel));
        if(node == NULL) return -1;
        s = strtok(ks->s, "\t");
        node->tid = bam_name2id(hdr, s);
        s = strtok(NULL, "\t");
        node->start = strtol(s, NULL, 10);
        s = strtok(NULL, "\n");
        node->end = strtoul(s, NULL, 10);
        insertNode(node);
    }
    if(ks->s) free(ks->s);
    free(ks);
    ks_destroy(kstr);
    return 0;
}

//Stub main
int main(int argc, char *argv[]) {
    htsFile *fp, *of;
    bam_hdr_t *hdr;
    gzFile bed;
    faidx_t *fai;

    //Open
    fp = sam_open(argv[1], "rb"); //Need to include an option to open a fasta file
    hdr = sam_hdr_read(fp);
    initTargetNodes();
    bed = gzopen(argv[2], "rb");
    assert(bed2list(bed, hdr) == 0); //This should be better handled
    gzclose(bed);
    fai = fai_load(argv[3]);

    //Write the output header
    of = sam_open(argv[4], "wb");
    sam_hdr_write(of, hdr);

    //Process the reads
    lastTargetNode = firstTargetNode->next;
    processReads(fp, hdr, of, 17, fai);

    //Clean up
    fprintf(stderr, "[main] 1\n"); fflush(stderr);
    destroyTargetNodes();
    fprintf(stderr, "[main] 2\n"); fflush(stderr);
    bam_hdr_destroy(hdr);
    fprintf(stderr, "[main] 3\n"); fflush(stderr);
    sam_close(fp);
    fprintf(stderr, "[main] 4\n"); fflush(stderr);
    sam_close(of);
    fprintf(stderr, "[main] 5\n"); fflush(stderr);
    fai_destroy(fai);
    fprintf(stderr, "[main] 6\n"); fflush(stderr);
    return 0;
}
