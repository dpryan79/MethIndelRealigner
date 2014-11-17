#include "realigner.h"
#include "htslib/faidx.h"
#include "SSW/ssw.h"

char int2base[32] = {0, 'A', 'T', 0, 'G', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N', \
                     0, 'A', 'C', 0, 'A', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N'};

                        // A CG  T  N
int8_t scoreMatrix[16] = { 2,-2,-2, 0,
                          -2, 2,-2, 0,
                          -2,-2, 2, 0,
                           0, 0, 0, 0};

//This is basically calculate_positions() from bison
//The output must be free()d
int32_t * bam2positions(bam1_t *b) {
    int32_t *positions, i = 0, j = 0;
    uint32_t *cigar = bam_get_cigar(b);
    int32_t prev = b->core.pos;
    int op, op_len, offset=0;

    positions = calloc(b->core.l_qseq, sizeof(int32_t));
    assert(positions);

    for(i=0; i<b->core.n_cigar; i++) {
        op=bam_cigar_op(cigar[i]);
        op_len = bam_cigar_oplen(cigar[i]);
        switch(op) {
        case 0 : //M
        case 7 : //=
        case 8 : //X
            for(j=0; j<op_len; j++) {
                assert(offset<b->core.l_qseq);
                positions[offset] = prev++;
                offset++;
            }
            break;
        case 1 : //I
        case 4 : //S
            //Assign the previous position
            for(j=0; j<op_len; j++) {
                assert(offset<b->core.l_qseq);
                positions[offset] = prev;
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
    int32_t *positions = NULL; //Could reuse this like a kstring_t
    int32_t refStart = (start-k>=0) ? start-k : 0;
    int32_t refEnd = refStart+refLen-1;
    uint64_t h;
#ifdef DEBUG
    char *tmp = calloc(k+1, sizeof(char));
    fprintf(stderr, "[bam2kmer] endpos is %"PRId32"\n", bam_endpos(b));
#endif

    positions = bam2positions(b);
    findPositions(positions, b->core.l_qseq, start, end-1, &start2, &end2);
    if(start2==-1) {
        free(positions);
        return;
    }
    //Does this even overlap the ROI?
    if(positions[start2]-refStart+end2-start2+1 <= 2*k) {
#ifdef DEBUG
        fprintf(stderr, "[bam2kmer] %s with bounds %"PRId32"-%"PRId32" doesn't overlap ROI from %"PRId32"-%"PRId32"\n", bam_get_qname(b), b->core.pos, bam_endpos(b), start,end);
#endif
        free(positions);
        return;
    }
#ifdef DEBUG
    fprintf(stderr, "[bam2kmer] start2 %"PRId32" end2 %"PRId32" refLen %"PRId32"\n", start2, end2, refLen);
    fprintf(stderr, "[bam2kmer] positions[start2] %"PRId32" start %"PRId32"\n", positions[start2], start);
    fprintf(stderr, "[bam2kmer] positions[end2] %"PRId32" end %"PRId32"\n", positions[end2], end);
#endif

    //Grow ks as needed
    if(ks->m < (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart))
        ks_resize(ks, (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart));
    ks->l = (refEnd-positions[end2])+(end2-start2+1)+(positions[start2]-refStart);

    //Add 5' reference sequence
#ifdef DEBUG
    fprintf(stderr, "[bam2kmer] copying first %"PRId32" chars onto 5' end\n", positions[start2]-refStart);
#endif
    if(positions[start2]-refStart) {
        if(offset==0) strncpy(ks->s, CT, positions[start2]-refStart);
        else strncpy(ks->s, GA, positions[start2]-refStart);
    }
#ifdef DEBUG
    ks->s[positions[start2]-refStart] = '\0';
    fprintf(stderr, "[bam2kmer] %s\n", ks->s);
    fprintf(stderr, "[bam2kmer] adding read from base %"PRId32" through %"PRId32" inclusive\n", start2, end2);
#endif

    //Extract the C->T and G->A sequences
    for(i=start2; i<=end2; i++) {
        ks->s[positions[start2]-refStart+i-start2] = int2base[offset+bam_seqi(bam_get_seq(b), i)];
    }
    ks->s[ks->l] = '\0';
#ifdef DEBUG
    ks->s[positions[start2]-refStart+i-start2] = '\0';
    fprintf(stderr, "[bam2kmer] %s\n", ks->s);
    fprintf(stderr, "[bam2kmer] current length should be %i\n", positions[start2]-refStart+i-start2);
    fprintf(stderr, "[bam2kmer] Final length should be %" PRId64"\n", ks->l);
#endif

    //Add 3' reference sequence
#ifdef DEBUG
    fprintf(stderr, "[bam2kmer] Copying over the last %"PRId32"\n", refEnd-positions[end2]);
    assert(refEnd-positions[end2] < refLen-k);
#endif
    if(refEnd-positions[end2]) {
        if(offset==0) strncpy(ks->s+positions[start2]-refStart+end2-start2+1,
                              CT+refLen-refEnd+positions[end2],
                              refEnd-positions[end2]);
        else strncpy(ks->s+positions[start2]-refStart+end2-start2+1,
                     GA+refLen-refEnd+positions[end2],
                     refEnd-positions[end2]);
    }
#ifdef DEBUG
    fprintf(stderr, "[bam2kmer] Final read sequence %s\n", ks->s); fflush(stderr);
#endif

    //Add the kmers, ignoring the first and last
    for(i=1; i<ks->l-k; i++) {
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
//Return -1 if there is no best alignment, which will result in the read not being realigned (we should mark these)
int32_t findBestAlignment(s_align **al, int32_t len, uint32_t *counts, int32_t readLen) {
    int32_t best=-1, bestScore = -1, i;
    uint32_t bestCount = 0;

    //Get the maximum score
    for(i=0; i<len; i++)
        if(al[i]->score1 > bestScore && al[i]->cigarLen == 1 && readLen == al[i]->read_end1-al[i]->read_begin1)
            bestScore = al[i]->score1;

    //Given a score, find the path with the highest count
    for(i=0; i<len; i++) {
        if(al[i]->score1 == bestScore && counts[i] > bestCount && al[i]->cigarLen == 1 && readLen == al[i]->read_end1-al[i]->read_begin1) {
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
            free(str->s);
            free(str);
        }
    }

    //Do we need to move everything after the cigar?
    if(b->core.n_cigar != len) {
        //Do we need to realloc b->data?
        if(b->m_data-b->l_data < 4*(len-b->core.n_cigar)) {
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
    int32_t oplen = 0, op=0,oplen2=0, op2=0, readPos = 0, newStartPos = -1;
    int32_t refPos = 0;
    int32_t i;

    assert(newCIGARArray);
    assert(readEndPos-readStartPos == refEndPos-refStartPos);

    /* We deal with things in 4 segments:
     * (1) read 5' overhang
     * (2) read 5' underhang
     * (3) overlap region
     * (4) read 3' overhang
     * Note that a read 3' underhang can be ignored!
     */
#ifdef DEBUG
    fprintf(stderr, "[updateAlignment] readStartPos %" PRId32 " readpos %" PRId32 " readEndPos %" PRId32 "\n", readStartPos, readPos, readEndPos);
    fprintf(stderr, "[updateAlignment] refStartPos %" PRId32 "\n", refStartPos);
#endif
    if(readStartPos > 0) { //5' overhang
        while(readPos < readStartPos) {
            if(oplen == 0) {
                assert(read_cigar_opnum < b->core.n_cigar);
                op = bam_cigar_op(bam_get_cigar(b)[read_cigar_opnum]);
                oplen = bam_cigar_oplen(bam_get_cigar(b)[read_cigar_opnum++]);
            }
            //If the op consumes the query and extends past the while-loop bounds, then trim
#ifdef DEBUG
            fprintf(stderr, "[updateAlignment] readPos %" PRId32 " oplen %"PRId32" readStartPos %"PRId32" bam_cigar_type(op) %i\n", readPos, oplen, readStartPos, bam_cigar_type(op));
#endif
            if(readPos + oplen>= readStartPos && (bam_cigar_type(op)&1)) {
                oplen = readStartPos-readPos;
                if(bam_cigar_type(op)&1) readPos += oplen;
                break;
            } else { //Append to newCIGARArray
                newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
                if(bam_cigar_type(op)&1) readPos += oplen;
                oplen = 0;
            }
        }
        //If refStartPos >0, we need to catch up to it
        if(refStartPos > 0) {
            while(refPos < refStartPos) {
                if(oplen2 == 0) {
                    assert(ref_cigar_opnum < al->cigarLen);
                    op2 = bam_cigar_op(al->cigar[ref_cigar_opnum]);
                    oplen2 = bam_cigar_oplen(al->cigar[ref_cigar_opnum++]);
                }
#ifdef DEBUG
                fprintf(stderr, "[updateAlignment] ref got cigar %"PRId32"%c\n", oplen2, BAM_CIGAR_STR[op2]);
#endif
                if(refPos + oplen2>= refStartPos && bam_cigar_type(op)&2) {
                    oplen2 = refPos+oplen2-refStartPos;
                    refPos = refStartPos;
                    break;
                }
                if(bam_cigar_type(op)&2) {
                    refPos += oplen2;
                    oplen2=0;
                }
            }
        }
#ifdef DEBUG
    fprintf(stderr, "[updateAlignment] 5' overhang CIGAR ");
    for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
    fprintf(stderr, "\n");
#endif
    } else if(refStartPos>0) { //5' Underhang, note new start position
        newStartPos = refStart;
#ifdef DEBUG
        fprintf(stderr, "[updateAlignment] refStart is %" PRId32 " so we start with newStartPos = that\n", refStart);
#endif
        while(refPos < refStartPos) {
            if(oplen2 == 0) {
                assert(ref_cigar_opnum < al->cigarLen);
                op2 = bam_cigar_op(al->cigar[ref_cigar_opnum]);
                oplen2 = bam_cigar_oplen(al->cigar[ref_cigar_opnum++]);
            }
#ifdef DEBUG
            fprintf(stderr, "[updateAlignment] refPos %" PRId32 " refStartPos %" PRId32 "\n", refPos, refStartPos);
            fprintf(stderr, "[updateAlignment] oplen %" PRId32 " op %c\n", oplen2, BAM_CIGAR_STR[op2]);
            fprintf(stderr, "[updateAlignment] bam_cigar_type(op) == %i\n", bam_cigar_type(op2));
            fflush(stderr);
#endif
            if(refPos + oplen2>= refStartPos && bam_cigar_type(op2)&1) {
                if(bam_cigar_type(op2) &2) newStartPos += refStartPos;
                oplen2 = refPos+oplen2-refStartPos;
#ifdef DEBUG
                fprintf(stderr, "[updateAlignment] Truncating oplen to %" PRId32 "\n", oplen2);
#endif
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
#ifdef DEBUG
        fprintf(stderr, "[updateAlignment] 5' underhang CIGAR ");
        for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
        fprintf(stderr, "\n");
#endif
    }
#ifdef DEBUG
    fprintf(stderr, "[updateAlignment] Remnant CIGAR is %" PRId32 "%c\n", oplen, BAM_CIGAR_STR[op]);
    fprintf(stderr, "[updateAlignment] Remnant CIGAR2 is %" PRId32 "%c\n", oplen2, BAM_CIGAR_STR[op2]);
    fflush(stderr);
    fprintf(stderr, "[updateAlignment] Before overlap region: refPos %" PRId32 " refEndPos %" PRId32 "\n", refPos, refEndPos);
#endif

    //Deal with the overlap region
    while(refPos <= refEndPos) {
        if(oplen2 == 0) {
            assert(ref_cigar_opnum < al->cigarLen);
            op2 = bam_cigar_op(al->cigar[ref_cigar_opnum]);
            oplen2 = bam_cigar_oplen(al->cigar[ref_cigar_opnum++]);
        }
#ifdef DEBUG
        fprintf(stderr, "[updateAlignment] 1 Found op %"PRId32"%c bam_cigar_type(op2) %i\n", oplen2, BAM_CIGAR_STR[op2], bam_cigar_type(op2));
        fprintf(stderr, "[updateAlignment] 1 refPos %" PRId32 " refEndPos+1 %" PRId32 "\n", refPos, refEndPos+1);
#endif
        if(refPos + oplen2>= refEndPos+1 && bam_cigar_type(op2)&1) {
            oplen2 = refEndPos-refPos+1;
#ifdef DEBUG
            fprintf(stderr, "[updateAlignment] Truncating the operation to length %" PRId32" and breaking\n", oplen2);
#endif
            if(oplen) {
                if(op == op2) {
                    oplen2 += oplen;
                } else newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
                oplen = 0;
            }
            refPos += oplen2;
            break;
        } else { //Append to newCIGARArray
            if(oplen) {
#ifdef DEBUG
                fprintf(stderr, "[updateAlignment] Comparing %c and %c\n", BAM_CIGAR_STR[op], BAM_CIGAR_STR[op2]);
#endif
                if(op == op2) {
#ifdef DEBUG
                    fprintf(stderr, "[updateAlignment] refPos before %" PRId32" ", refPos);
#endif
                    if(bam_cigar_type(op2) &1) refPos += oplen2;
#ifdef DEBUG
                    fprintf(stderr, " and %" PRId32" after\n", refPos);
#endif
                    oplen2 += oplen;
                } else {
                    newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
                    if(bam_cigar_type(op2) &1) refPos += oplen2;
                }
                oplen = 0;
            } else {
                if(bam_cigar_type(op2) &1) refPos += oplen2;
            }
#ifdef DEBUG
            fprintf(stderr, "[updateAlignment] pushing %"PRId32"%c\n", oplen2, BAM_CIGAR_STR[op2]);
#endif
            newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op2, oplen2);
            oplen2 = 0;
        }
        oplen2 = 0;
#ifdef DEBUG
        fprintf(stderr, "[updateAlignment] 2 refPos %" PRId32 " refEndPos+1 %" PRId32 "\n", refPos, refEndPos+1);
#endif
    }
#ifdef DEBUG
    fprintf(stderr, "[updateAlignment] Remnant CIGAR is %" PRId32 "%c\n", oplen2, BAM_CIGAR_STR[op2]);
#endif
    oplen = 0;
#ifdef DEBUG
    fprintf(stderr, "[updateAlignment] Post-overlap CIGAR ");
    for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
    fprintf(stderr, "\n");
    fprintf(stderr, "[updateAlignment] Before 3' overhang, readPos %" PRId32 "readEndPos %" PRId32 " b->core.l_qseq-1 %" PRId32"\n", readPos, readEndPos, b->core.l_qseq-1);
#endif

    //3' Overhang
    if(readEndPos < b->core.l_qseq-1) {
        //Iterate up to the relevant CIGAR operation
        oplen = 0;
        readPos = -1;
        for(read_cigar_opnum=0; read_cigar_opnum<b->core.n_cigar; read_cigar_opnum++) {
            op = bam_cigar_op(bam_get_cigar(b)[read_cigar_opnum]);
            oplen = bam_cigar_oplen(bam_get_cigar(b)[read_cigar_opnum]);
#ifdef DEBUG
            fprintf(stderr, "[updateAlignment] Found %"PRId32"%c\n", oplen, BAM_CIGAR_STR[op]);
#endif
            if(bam_cigar_type(op) & 1) {
                if(readPos+oplen>readEndPos) {
                    oplen = oplen+readPos-readEndPos;
#ifdef DEBUG
                    fprintf(stderr, "[updateAlignment] Truncating to %"PRId32" due to readPos=%"PRId32" and readEndPos %"PRId32"\n", oplen, readPos, readEndPos);
#endif
                    readPos += oplen;
                    read_cigar_opnum++; //Otherwise, the break will prevent this
                    break;
                }
                readPos += oplen;
            }
        }
#ifdef DEBUG
        fprintf(stderr, "[updateAlignment] read_cigar_opnum %"PRId32"\n", read_cigar_opnum);
        fprintf(stderr, "[updateAlignment] remnant CIGAR op is %"PRId32"%c\n", oplen, BAM_CIGAR_STR[op]);
        fprintf(stderr, "[updateAlignment] remnant CIGAR2 op is %"PRId32"%c\n", oplen2, BAM_CIGAR_STR[op2]);
        fprintf(stderr, "[updateAlignment] 3' overhang remnant CIGAR ");
        for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
        fprintf(stderr, "\n");
#endif
        //Can we merge this OP with that from the overlap?
        if(oplen2) {
            if(op==op2) {
                oplen +=oplen2;
            } else {
                newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op2, oplen2);
                oplen2 = 0;
            }
        }
#ifdef DEBUG
        fprintf(stderr, "[updateAlignment] remnant CIGAR op is %"PRId32"%c\n", oplen, BAM_CIGAR_STR[op]);
        fprintf(stderr, "[updateAlignment] 3' overhang CIGAR ");
        for(i=0; i<n_cigar; i++) fprintf(stderr, "%"PRId32"%c", bam_cigar_oplen(newCIGARArray[i]), BAM_CIGAR_STR[bam_cigar_op(newCIGARArray[i])]);
        fprintf(stderr, "\n");
#endif

        //Push the remnant CIGAR, if it exists
        if(oplen) {
#ifdef DEBUG
            fprintf(stderr, "[updateAlignment] Pushing remnant %"PRId32"%c\n", oplen, BAM_CIGAR_STR[op]);
#endif
            newCIGARArray = pushCIGAR(newCIGARArray, n_cigar++, &max_cigar, op, oplen);
        }
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
    assert(n_cigar>0);
    //Ensure that the new CIGAR string matchs b->core.l_qseq
    oplen2 = 0;
    for(i=0; i<n_cigar;i++) {
        op = bam_cigar_op(newCIGARArray[i]);
        oplen = bam_cigar_oplen(newCIGARArray[i]);
        if(bam_cigar_type(op)&1) oplen2+=oplen;
#ifdef DEBUG
        fprintf(stderr, "[updateAlignment] Found %"PRId32"%c, type %i oplen2 %"PRId32"\n", oplen, BAM_CIGAR_STR[op], bam_cigar_type(op), oplen2);
    }
    fflush(stderr);
#else
    }
#endif
    assert(oplen2==b->core.l_qseq);
    updateCigar(b, newCIGARArray, n_cigar);
    free(newCIGARArray);
    return b;
}

//Align all of the paths to a reference sequence
s_align **alignPaths2Ref(paths *p, int8_t *ref, int32_t refLen, char *refSeq, int32_t *refIndex, int32_t k) {
    int32_t i;
    s_align **alignments = malloc(sizeof(s_align*) * p->l);
    assert(alignments);
    for(i=0; i<p->l; i++) {
        alignments[i] = SemiGlobalAlignment(ref, refLen, p->conv[i], p->len[i], k);
        if(*refIndex == -1 && alignments[i]->cigarLen == 1)
            if(strcmp(refSeq, p->path[i]) == 0) *refIndex = i;
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
    findPositions(positions, b->core.l_qseq, refLBound, refRBound-1, readLBound, readRBound);
    if(*readLBound == -1) { //Read doesn't overlap ROI
        free(readal);
        free(positions);
        return NULL;
    }
    subreadL = (*readRBound)-(*readLBound)+1;
    if(subreadL >= *subreadM) { //Grow subreadSeq if needed
        (*subreadM) = subreadL+1;
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
        fprintf(stderr, "[alignReads2Paths] readal[%i]->ref_begin1 %" PRId32 " ref_end1 %" PRId32 "\n", i, readal[i]->ref_begin1, readal[i]->ref_end1);
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
    CTal = alignPaths2Ref(CTpaths, CTref, refLen, CT, &(refIndex[0]), k);
    GAal = alignPaths2Ref(GApaths, GAref, refLen, GA, &(refIndex[1]), k);

    //Iterate over the reads, aligning them to the paths
    refLBound = ((*heap)->start-k >= 0) ? (*heap)->start-k : 0;
    refRBound = refLBound + refLen-2*k;
    for(i=0; i<(*heap)->l; i++) {
        strand = getStrand((*heap)->heap[i]);
#ifdef DEBUG
        fprintf(stderr, "[realignHeapCore] processing heap->heap[%i]\n", i);
#endif
        readal[i] = alignReads2Paths((*heap)->heap[i], strand, &subreadM, &subreadSeq, (strand==1) ? CTpaths : GApaths, refLBound, refRBound, readLBound+i, readRBound+i);
    }

    //Count number of best/equally best alignments/path
    uint32_t *CTcounts = calloc(CTpaths->l, sizeof(uint32_t));
    uint32_t *GAcounts = calloc(GApaths->l, sizeof(uint32_t));
    for(i=0; i<(*heap)->l; i++) {
        strand = getStrand((*heap)->heap[i]);
        if(readal[i] == NULL) continue;
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

    //Pick a best path for each alignment given all of the others and update it
    refLBound = ((*heap)->start-2*k >= 0) ? (*heap)->start-2*k : 0; //The previous bounds were only for extracting subreads
    for(i=0; i<(*heap)->l; i++) {
        strand = getStrand((*heap)->heap[i]);
        if(readal[i] == NULL) continue;
        if(strand&1) {
            bestAl = findBestAlignment(readal[i], CTpaths->l, CTcounts, readRBound[i]-readLBound[i]);
        } else {
            bestAl = findBestAlignment(readal[i], GApaths->l, GAcounts, readRBound[i]-readLBound[i]);
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
        else if(!(strand&1) && refIndex[1] != bestAl)
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
        if(readal[i] == NULL) continue;
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
    free(CTal);
    free(GAal);
    free(subreadSeq);
    free(readLBound);
    free(readRBound);
}

paths *addRefPath(char *Seq, int len, paths *p) {
    int32_t i, found = 0;
    for(i=0; i<p->l; i++) {
        if(p->len[i] == len) {
            if(strcmp(p->path[i], Seq) == 0) {
                found = 1;
                break;
            }
        }
    }
    if(found==0) {
#ifdef DEBUG
        fprintf(stderr, "[addRefPath] Manually added ref path %s\n", Seq); fflush(stderr);
#endif
        p->l++;
        if(p->l >= p->m) {
            p->m = p->l;
            p->len = realloc(p->len, sizeof(int32_t)*(p->m));
            assert(p->len);
            p->path = realloc(p->path, sizeof(char *)*(p->m));
            assert(p->path);
        }
        p->len[p->l-1] = len;
        p->path[p->l-1] = strdup(Seq);
    }
    return p;
}

//k is the kmer
void realignHeap(alignmentHeap *heap, int k, faidx_t *fai) {
    bf *filter = bf_init(heap->end-heap->start, k);
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
    fprintf(stderr, "[realignHeap] heap->start %"PRId32" heap->end %"PRId32"\n", heap->start, heap->end);
    fprintf(stderr, "[realignHeap] fetching sequence from %"PRId32"-%"PRId32" of length %i\n", start2, end+k-1,len);
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

    //Ensure that the reference path is included (it won't be if it contains a cycle!)
    CTpaths = addRefPath(CT, len, CTpaths);
    GApaths = addRefPath(GA, len, GApaths);

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
