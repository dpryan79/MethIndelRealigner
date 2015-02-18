#include "realigner.h"
#include "htslib/kseq.h"
#include "htslib/faidx.h"
#include "htslib/bgzf.h"
#include "threads.h"
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

void processReads(htsFile *fp, bam_hdr_t *hdr, htsFile *of, int k, faidx_t *fai, int depth, int nt, int quiet, int threshold) {
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
            if(!quiet) {
                fprintf(stderr, "ROI %s:%"PRId32"-%"PRId32"\n", GLOBAL_HEADER->target_name[bounds->tid], bounds->start, bounds->end); fflush(stderr);
            }
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
                fprintf(stderr, "[processReads] b->core.pos %" PRId32 " endpos %" PRId32 "\n", b->core.pos, bam_endpos(b)-1);
                fprintf(stderr, "[processReads] start %"PRId32" end %"PRId32" b->core.pos %" PRId32" %s\n", heap->start, heap->end, b->core.pos, bam_get_qname(b));
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
                if(heap->l) realignHeap(heap, k, fai, nt, threshold, quiet);
#ifdef DEBUG
                fprintf(stderr, "[processReads] Last heap realigned\n"); fflush(stderr);
#endif
                writeHeap(of, hdr, heap, 1);
                heap->l = 0;
                break;
            } else if(status == 2) {
                if(quiet < 2) fprintf(stderr, "[processReads] Skipping %s:%"PRId32"-%"PRId32", too many reads\n", hdr->target_name[heap->heap[0]->core.tid], heap->start, heap->end);
            } else {
#ifdef DEBUG
                fprintf(stderr, "[processReads] %s is outside the ROI, realigning the heap\n", bam_get_qname(b)); fflush(stderr);
#endif
                realignHeap(heap, k, fai, nt, threshold, quiet);
            }
            lastTargetNode = lastTargetNode->next; //Move to the next ROI
            heap = writeHeapUntil(of, hdr, heap, b, fp);
            if(heap->l) {//The heap wasn't flushed
                if(!quiet) {
                    fprintf(stderr, "ROI %s:%"PRId32"-%"PRId32"\n", GLOBAL_HEADER->target_name[lastTargetNode->tid], lastTargetNode->start, lastTargetNode->end); fflush(stderr);
                }
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
int bed2list(gzFile fp, bam_hdr_t *hdr, int32_t k, int quiet) {
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
        if(lastTargetNode != firstTargetNode) {
            if(lastTargetNode->tid == node->tid && lastTargetNode->end+2*k > node->start) {
                lastTargetNode->end = node->end;
                free(node);
                continue;
            }
        }
        insertAfter(node, lastTargetNode);
        lastTargetNode = lastTargetNode->next;
        total++;
    }
    lastTargetNode = firstTargetNode;
    if(ks->s) free(ks->s);
    free(ks);
    ks_destroy(kstr);
    if(!quiet) {
        fprintf(stderr, "Found %"PRIu32" ROIs\n", total); fflush(stderr);
    }
    return 0;
}

void realigner_usage() {
    fprintf(stderr, "Usage: MethIndelRealigner realign [OPTIONS] <input.bam> <reference.fa> [output.bam]\n");
    fprintf(stderr, 
"\nNote that input.bam must be coordinate-sorted and indexed. Similarly,\n"
"reference.fa must be indexed with samtools faidx. If an output file name isn't\n"
"given, then the output will be written to stdout (so use '>' to redirect to a\n"
"file).\n"
"\nThe input and output formats will always be the same. So BAM input produces\n"
"BAM output and CRAM input produces CRAM output.\n"
"\nOPTIONS\n"
"-C       Produce results in CRAM format, regardless of whether the input is in\n"
"         CRAM format or not.\n"
"-q INT   The minimum MAPQ value to process an alignment in the target creation\n"
"         step or to realign it in the realignment step. The default is 10.\n"
"-k INT   The k-mer size to use. If you manually ran the TargetCreator, the\n"
"         setting you used there should match. This must be an odd number.\n"
"         Default is 25.\n"
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
"--quiet  Suppress printing the ROI (region of interest) that's being processed.\n"
"         Specifying this twice will suppress everything except error messages.\n"
"--index  If the output is being written directly to a file (rather than to the\n"
"         screen), attempt to index it upon completion.\n"
"\n"
"Advanced options:\n"
"         Unless you're familiar with how graph algorithms work, it's unwise to\n"
"         change these.\n"
"\n"
"--nt INT Number of threads to use for realignment of graph paths to the\n"
"         reference sequence and reads to the paths. While increasing the\n"
"         number of threads will speed up processing, due to how things are\n"
"         implemented there will also be a lot of wasted computer resources (the\n"
"         threads use spin locks because context switching overhead dwarfs thread\n"
"         execution time). Particularly on shared systems, this may not be ideal.\n"
"         Also note that while using --nt 1 will request only a single thread,\n"
"         using --nt 2 will produce 4 independent processing threads that are\n"
"         always running (so, a total of 5 concurrent threads, and --nt 3 would\n"
"         result in 7 threads). The default value is 1.\n"
"--minKmerCount INT The minimum number of times a Kmer must appear for it to be\n"
"         used in finding paths through the de Bruijn graph. Lower numbers will\n"
"         produce more (possibly too many, see the --breadth setting) paths,\n"
"         though a given alignment is more likely to then align perfectly to it.\n"
"         The default is 4 and this value must be >0 and not more than 15.\n"
"--maxInsert INT The maximum insert size that will be looked for. The default is\n"
"         0, meaning no limit is used. Setting this to your maximum read length\n"
"         can be useful in repeat regions where there are many misaligned reads.\n"
"--breadth INT Maximum breadth of the C->T or G->A graph. If more paths are\n"
"         found in either graph than this, then no realignment will take place.\n"
"         This effectively limits realignment around very low complexity regions\n"
"         where the amount of memory and time needed grows exponentially. Setting\n"
"         this to 0 sets no maximum. This value must be >= 0. The default is 300.\n"
"--noCycles If you specify this flag, then cycles in paths traversing an InDel\n"
"         will be avoided. While this is a very good idea in theory (it speeds\n"
"         things up by preventing alignments to aberrant repeat sequences), in\n"
"         practice the low information content of bisulfite-converted DNA makes\n"
"         the presence of cycles very common. While the a path representing the\n"
"         reference sequence will always be present (regardless of it having\n"
"         cycles), if you specify --noCycles and happen to choose a setting for\n"
"         -k that's too low to accurately represent a given path without cycles\n"
"         then it won't be traversed. While current alignments to that path will\n"
"         likely be maintained (since realignment will increase the edit\n"
"         distance, producing realignments that will be rejected in favor of the\n"
"         initial alignments), new alignments will be impossible. Note that the\n"
"         maximum and minimum path length is dictated by the alignments\n"
"         traversing an ROI, so and endless number of cycles won't be followed\n"
"         when creating paths in any case.\n");
}

int realigner_main(int argc, char *argv[]) {
    htsFile *fp, *of;
    bam_hdr_t *hdr;
    gzFile bed;
    faidx_t *fai;
    int c, bedSet=0, depth = 1000, kmer = 25, quiet = 0;
    int ROIdepth = 4, nCompThreads = 1, nt = 1, threshold = 4;
    int doIndex = 0, CRAM = 0;
    uint32_t total = 0;
    MAXBREADTH = 300;
    MAXINSERT = 0;
    MINMAPQ = 10;
    NOCYCLES = 0;

    static struct option lopts[] = {
        {"help",        0, NULL,'h'},
        {"ROIdepth",    1, NULL,'d'},
        {"breadth",     1, NULL, 1},
        {"maxInsert",   1, NULL, 2},
        {"noCycles",    0, NULL, 3},
        {"nt",          1, NULL, 4},
        {"quiet",       0, NULL, 5},
        {"minKmerCount",1, NULL, 6},
        {"index",       0, NULL, 7},
        {NULL,          0, NULL, 0}
    };

    opterr = 0; //Disable error messages
    while((c = getopt_long(argc, argv, "Ck:l:D:@:q:d:h", lopts, NULL)) >= 0) {
        switch(c) {
        case 'C' :
            CRAM = 1;
            break;
        case 'h' :
            realigner_usage();
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
            MINMAPQ = atoi(optarg);
            break;
        case '@' :
            nCompThreads = atoi(optarg);
            break;
        case 1 :
            MAXBREADTH = atoi(optarg);
            if(MAXBREADTH<0) {
                fprintf(stderr, "--breadth must be >=0! Resetting it to 0.\n");
                MAXBREADTH = 0;
            }
            break;
        case 2 :
            MAXINSERT = atoi(optarg);
            if(MAXINSERT<0) {
                fprintf(stderr, "--maxInsert must be >=0! Resetting it to 0.\n");
                MAXINSERT = 0;
            }
            break;
        case 3 :
            NOCYCLES = 1;
            break;
        case 4 :
            nt = atoi(optarg);
            if(nt<1) nt = 1;
            break;
        case 5 :
            quiet += 1;
            break;
        case 6 :
            threshold = atoi(optarg);
            if(threshold<1) threshold = 1;
            if(threshold>15) threshold = 15;
            break;
        case 7:
            doIndex = 1;
            break;
        case '?' :
        default :
            fprintf(stderr, "Invalid option '%s'\n", argv[optind-1]);
            realigner_usage();
            return 1;
            break;
        }
    }

    if(argc == 1) {
        fprintf(stderr, "A\n");
        realigner_usage();
        return 0;
    }
    if(argc-optind != 3 && argc-optind != 2) {
        fprintf(stderr, "B\n");
        realigner_usage();
        return 1;
    }

    //Open the input files
    fp = sam_open(argv[optind], "r");
    hts_set_fai_filename(fp, argv[optind+1]);
    hdr = sam_hdr_read(fp);
    initTargetNodes();

    if(bedSet) {
        assert(bed2list(bed, hdr, kmer, quiet) == 0); //This should be better handled
        gzclose(bed);
    } else {
        //Run TargetCreator ourselves
        initTargetNodes();
        findInDels(fp, hdr, MINMAPQ, 5);
        total = depthFilter(ROIdepth);
        if(!quiet) {
            fprintf(stderr, "Found %"PRIu32" ROIs\n", total);
            fflush(stderr);
        }
        bam_hdr_destroy(hdr);
        sam_close(fp);
        //Reopen the input BAM file
        fp = sam_open(argv[optind], "rb");
        hts_set_fai_filename(fp, argv[optind+1]);
        hdr = sam_hdr_read(fp);
    }

    fai = fai_load(argv[optind+1]);
    GLOBAL_FAI = fai;
    if(argc-optind==2) {
        //write to stdout
        if(fp->is_cram || CRAM) of = sam_open("-", "wc");
        else of = sam_open("-", "wb0");
    } else {
        if(fp->is_cram || CRAM) of = sam_open(argv[optind+2], "wc");
        else of = sam_open(argv[optind+2], "wb");
    }
    hts_set_fai_filename(fp, argv[optind+1]);
    hts_set_fai_filename(of, argv[optind+1]);
    if(nCompThreads>1) hts_set_threads(of, nCompThreads);
    sam_hdr_write(of, hdr);
    GLOBAL_HEADER = hdr;

    //setup the condition variables
    if(nt>1) setup_threads(nt);

    //Process the reads
    lastTargetNode = firstTargetNode->next;
    processReads(fp, hdr, of, kmer, fai, depth, nt, quiet, threshold);

    //Clean up
    destroyTargetNodes();
    bam_hdr_destroy(hdr);
    sam_close(fp);
    sam_close(of);
    fai_destroy(fai);
    destroy_threads(nt);

    //If the output file was not stdout, then try to index it
    if(argc-optind==3 && doIndex) {
        if(bam_index_build(argv[optind+2], 0) != 0) {
            fprintf(stderr, "Couldn't index %s, please manually sort and index it with samtools.\n", argv[optind+2]);
        }
    }

    return 0;
}
