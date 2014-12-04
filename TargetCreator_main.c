#include "realigner.h"
#include <getopt.h>

void usage(char *prog) {
    fprintf(stderr, "%s [OPTIONS] file.bam [output.bed]\n", prog);
    fprintf(stderr,
"\nNote that regions of interest (ROIs) within 5 bases of each other will be\n"
"merged, since they likely arise from the same event.\n"
"\nIf output.bed isn't specified, it'll be printed to stdout (the screen)\n"
"\nOPTIONS:\n"
"-q INT   The minimum MAPQ value to process an alignment in the target creation\n"
"         step or to realign it in the realignment step. The default is 10.\n"
"--ROIdepth INT Minimum depth covering a putative ROI to include it. This option\n"
"         is ignored if you supply a BED file. The default is 4, meaning that\n"
"         you need 4 reads supporting an InDel for realignment to occur around\n"
"         it.\n");
}

int main(int argc, char *argv[]) {
    htsFile *fp = NULL;
    bam_hdr_t *hdr;
    uint32_t total = 0;
    FILE *of = NULL;
    int ROIdepth = 4, c;
    MINMAPQ = 10;

    static struct option lopts[] = {
        {"help",     0, NULL, 'h'},
        {"ROIdepth", 1, NULL, 'd'}};

    while((c = getopt_long(argc, argv, "q:h", lopts, NULL)) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
            break;
        case 'd' :
            ROIdepth = atoi(optarg);
            break;
        case 'q' :
            MINMAPQ = atoi(optarg);
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
    if(argc-optind != 2 || argc-optind != 1) {
        usage(argv[0]);
        return 1;
    }

    //Open
    fp = sam_open(argv[optind], "rb");
    if(!fp) return 1;
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
