#include "realigner.h"
#include <getopt.h>

void TargetCreator_usage() {
    fprintf(stderr, "MethIndelRealigner TargetCreator [OPTIONS] <file.bam> [output.bed]\n");
    fprintf(stderr,
"\nNote that regions of interest (ROIs) within 5 bases of each other will be\n"
"merged, since they likely arise from the same event.\n"
"\nIf output.bed isn't specified, it'll be printed to stdout (the screen)\n"
"\nOPTIONS:\n"
"-f FILE  If the input is in CRAM format, specify where the reference sequence\n"
"         fasta file is located.\n"
"-q INT   The minimum MAPQ value to process an alignment in the target creation\n"
"         step or to realign it in the realignment step. The default is 10.\n"
"--ROIdepth INT Minimum depth covering a putative ROI to include it. This option\n"
"         is ignored if you supply a BED file. The default is 4, meaning that\n"
"         you need 4 reads supporting an InDel for realignment to occur around\n"
"         it.\n");
}

int TargetCreator_main(int argc, char *argv[]) {
    htsFile *fp = NULL;
    bam_hdr_t *hdr;
    char *ref = NULL;
    uint32_t total = 0;
    FILE *of = NULL;
    int ROIdepth = 4, c;
    MINMAPQ = 10;

    static struct option lopts[] = {
        {"help",     0, NULL, 'h'},
        {"ROIdepth", 1, NULL, 'd'},
        {NULL,       0, NULL,  0 }
    };

    opterr = 0;
    while((c = getopt_long(argc, argv, "q:f:h", lopts, NULL)) >= 0) {
        switch(c) {
        case 'h' :
            TargetCreator_usage();
            return 0;
            break;
        case 'd' :
            ROIdepth = atoi(optarg);
            break;
        case 'q' :
            MINMAPQ = atoi(optarg);
            break;
        case 'f' :
            ref = optarg;
            break;
        default :
            fprintf(stderr, "Invalid option '%s'\n", argv[optind-1]);
            TargetCreator_usage();
            return 1;
        }
    }

    if(argc == 1) {
        TargetCreator_usage();
        return 0;
    }
    if(argc-optind != 2 && argc-optind != 1) {
        TargetCreator_usage();
        return 1;
    }

    //Open
    fp = sam_open(argv[optind], "rb");
    if(!fp) return 1;
    if(ref) hts_set_fai_filename(fp, ref);
    if(argc-optind == 2) {
        of = fopen(argv[optind+1], "w");
        if(!of) return 1;
    } else {
        of = stdout;
    }
    hdr = sam_hdr_read(fp);
    if(!hdr) return 1;
    initTargetNodes();
    //Iterate
    findInDels(fp, hdr, MINMAPQ, 5);
    //Filter
    total = depthFilter(ROIdepth);
    //Write output
    writeTargets(of, hdr);
    //Clean up
    destroyTargetNodes();
    bam_hdr_destroy(hdr);
    sam_close(fp);
    if(argc-optind == 2) fclose(of);

    fprintf(stderr, "Found %" PRIu32 " sights.\n", total);
    return 0;
}
