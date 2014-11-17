#include "realigner.h"
#include "htslib/kseq.h"
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
#include <zlib.h>
#include <getopt.h>

KSTREAM_INIT(gzFile, gzread, 4096)

bam_hdr_t *GLOBAL_HEADER; //This just makes printing convenient

struct bounds {
    int32_t tid, start, end;
};

int regionNodeCmp(int32_t tid, int32_t start, int32_t end) {
    if(tid == lastTargetNode->tid) {
        if((start <= lastTargetNode->start && end >= lastTargetNode->start) || \
            (start >= lastTargetNode->start && start <= lastTargetNode->end)) return 0;
        else if(end < lastTargetNode->start) return end-lastTargetNode->start;
        assert(start > lastTargetNode->end);
        return start - lastTargetNode->end;
    } else return tid - lastTargetNode->tid;
}

int overlapsRegionsCore(int32_t tid, int32_t start, int32_t end) {
    int direction;

    if(lastTargetNode==NULL) return 0;

    direction = regionNodeCmp(tid,start,end);
    if(direction == 0) return 1;
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

void processReads(htsFile *fp, bam_hdr_t *hdr, htsFile *of, int k, faidx_t *fai, int depth) {
    bam1_t *b = bam_init1();
    struct bounds *bounds;
    alignmentHeap *heap = alignmentHeap_init(depth);
    int status; //1: EOF, 2: heap max, 3: Past ROI

    while(sam_read1(fp, hdr, b) > 0) {
#ifdef DEBUG
        fprintf(stderr, "[processReads] Found %s\n", bam_get_qname(b)); fflush(stderr);
#endif
        if(!(b->core.flag&4) && (bounds = overlapsRegion(b)) != NULL) {
#ifdef DEBUG
            fprintf(stderr, "[processReads] %s in region\n", bam_get_qname(b)); fflush(stderr);
#endif
            fprintf(stderr, "%s:%"PRId32"-%"PRId32"\n", GLOBAL_HEADER->target_name[bounds->tid], bounds->start, bounds->end); fflush(stderr);
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
            while(1) {
                if(sam_read1(fp, hdr, b)<0) {
                    status = 1;
                    break;
                }
#ifdef DEBUG
                fprintf(stderr, "[processReads] Found %s while in heap of length %i\n", bam_get_qname(b), heap->l); fflush(stderr);
                fprintf(stderr, "b->core.pos %" PRId32 " endpos %" PRId32 "\n", b->core.pos, bam_endpos(b)-1);
                fprintf(stderr, "start %"PRId32" end %"PRId32" b->core.pos %" PRId32" %s\n", heap->start, heap->end, b->core.pos, bam_get_qname(b));
#endif
inheap:         if(b->core.pos < heap->end && b->core.tid == heap->heap[0]->core.tid) {
                    //If the alignment and the first in the heap have the same pos...
                    if(b->core.pos == heap->heap[0]->core.pos) {
                        if((bounds = overlapsRegion(b)) == NULL) {
                            sam_write1(of, hdr, b);
                            continue;
                        }
                        free(bounds);
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
            if(status==2){ fprintf(stderr, "Before heap->l %" PRId32"\n", heap->l); fflush(stderr);}
            heap = writeHeapUntil(of, hdr, heap, depth);
            if(status==2){ fprintf(stderr, "After heap->l %" PRId32"\n", heap->l); fflush(stderr);}
            if(heap->l) {//The heap wasn't flushed
                fprintf(stderr, "%s:%"PRId32"-%"PRId32"\n", GLOBAL_HEADER->target_name[lastTargetNode->tid], lastTargetNode->start, lastTargetNode->end); fflush(stderr);
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
    uint32_t total = 0;

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
        insertAfter(node, lastTargetNode);
        lastTargetNode = lastTargetNode->next;
        total++;
    }
    lastTargetNode = firstTargetNode;
    if(ks->s) free(ks->s);
    free(ks);
    ks_destroy(kstr);
    fprintf(stderr, "Found %"PRIu32" ROIs\n", total); fflush(stderr);
    return 0;
}

void usage(char *prog) {
    fprintf(stderr, "Usage: %s [OPTIONS] input.bam reference.fa output.bam\n", prog);
    fprintf(stderr, 
"\nNote that input.bam must be coordinate-sorted and indexed. Similarly,\n"
"reference.fa must be indexed with samtools faidx.\n"
"\nOPTIONS\n"
"-k INT   The k-mer size to use. If you manually ran the TargetCreator, the\n"
"         setting you used there should match. This must be an odd number.\n"
"         Default is 17.\n"
"-l FILE  A (possibly gzipped) BED file listing regions to realign around. A\n"
"         region +/- k will be realigned. It is suggested that you run the\n"
"         TargetCreator and supply the resulting BED file here.\n"
"-D INT   Maximum heap depth. This effectively limits the depth covering any\n"
"         position that needs to be realigned. The default is 1000.\n"
"-@ INT   Number of compression threads (equivalent to the -@ option in\n"
"         samtools). Default 1.\n"
"--ROIdepth INT Minimum depth covering a putative ROI to include it. This option\n"
"         is ignored if you supply a BED file. The default is 4, meaning that\n"
"         you need 4 reads supporting an InDel for realignment to occur around\n"
"         it.\n"
"--ROIqual INT Minimum MAPQ score needed to process an alignment in the target\n"
"         creation step. This must be >0. The default is 10.\n");
}

//This should be exanded upon to attempt to index the fasta and BAM files if needed
int main(int argc, char *argv[]) {
    htsFile *fp, *of;
    bam_hdr_t *hdr;
    gzFile bed;
    faidx_t *fai;
    int c, bedSet=0, depth = 1000, kmer = 17;
    int ROIdepth = 4, ROIqual = 10, nCompThreads = 1;
    uint32_t total = 0;

    static struct option lopts[] = {
        {"help",     0, NULL, 'h'},
        {"ROIdepth", 1, NULL, 'd'},
        {"ROIqual",  1, NULL, 'q'}
    };

    while((c = getopt_long(argc, argv, "k:l:D:@:", lopts, NULL)) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
            break;
        case 'k' :
            kmer = atoi(optarg);
            if(!(kmer&1)) {
                fprintf(stderr, "-k must be an odd number, changing to %i\n", ++kmer);
            }
            break;
        case 'l' :
            bedSet=1;
            bed = gzopen(optarg, "rb");
            break;
        case 'D' :
            depth = atoi(optarg);
            break;
        case 'd' :
            ROIdepth = atoi(optarg);
            break;
        case 'q' :
            ROIqual = atoi(optarg);
            break;
        case '@' :
            nCompThreads = atoi(optarg);
            break;
        default :
            fprintf(stderr, "Invalid option '%c'\n", c);
            usage(argv[0]);
            return 1;
        }
    }

    if(argc == 1) {
        usage(argv[0]);
        return 0;
    }
    if(argc-optind != 3) {
        usage(argv[0]);
        return 1;
    }

    //Open the input files
    fp = sam_open(argv[optind], "rb");
    hdr = sam_hdr_read(fp);
    initTargetNodes();

    if(bedSet) {
        assert(bed2list(bed, hdr) == 0); //This should be better handled
        gzclose(bed);
    } else {
        //Run TargetCreator ourselves
        initTargetNodes();
        findInDels(fp, hdr, ROIqual, kmer);
        total = depthFilter(ROIdepth);
        fprintf(stderr, "Found %"PRIu32" ROIs\n", total+1);
        fflush(stderr);
        bam_hdr_destroy(hdr);
        sam_close(fp);
        //Reopen the input BAM file
        fp = sam_open(argv[optind], "rb");
        hdr = sam_hdr_read(fp);
    }

    fai = fai_load(argv[optind+1]);
    of = sam_open(argv[optind+2], "wb");
    if(nCompThreads>1) bgzf_mt(of->fp.bgzf, nCompThreads, 256);
    sam_hdr_write(of, hdr);
    GLOBAL_HEADER = hdr;

    //Process the reads
    lastTargetNode = firstTargetNode->next;
    processReads(fp, hdr, of, kmer, fai, depth);

    //Clean up
    destroyTargetNodes();
    bam_hdr_destroy(hdr);
    sam_close(fp);
    sam_close(of);
    fai_destroy(fai);
    return 0;
}

