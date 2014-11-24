#include "realigner.h"
#include <getopt.h>

void usage(char *prog) {
    fprintf(stderr, "%s [OPTIONS] file.bam output.bed\n", prog);
    fprintf(stderr,
"\nOPTIONS:\n"
"-q INT   The minimum MAPQ value to process an alignment in the target creation\n"
"         step or to realign it in the realignment step. The default is 10.\n"
"-k INT   The k-mer size to use. This must be an odd number. The default is 17.\n"
"--ROIdepth INT Minimum depth covering a putative ROI to include it. This option\n"
"         is ignored if you supply a BED file. The default is 4, meaning that\n"
"         you need 4 reads supporting an InDel for realignment to occur around\n"
"         it.\n");
}

int main(int argc, char *argv[]) {
    htsFile *fp;
    bam_hdr_t *hdr;
    uint32_t total = 0;
    FILE *of;
    int ROIdepth = 4, kmer = 17, c;
    MINMAPQ = 10;

    static struct option lopts[] = {
        {"help",     0, NULL, 'h'},
        {"ROIdepth", 1, NULL, 'd'}};

    while((c = getopt_long(argc, argv, "k:q:h", lopts, NULL)) >= 0) {
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
    if(argc-optind != 2) {
        usage(argv[0]);
        return 1;
    }

    //Open
    fp = sam_open(argv[optind], "rb");
    of = fopen(argv[optind+1], "w");
    hdr = sam_hdr_read(fp);
    initTargetNodes();
    //Iterate
    findInDels(fp, hdr, MINMAPQ, kmer);
    //Filter
    total = depthFilter(ROIdepth);
    //Write output
    writeTargets(of, hdr);
    //Clean up
    destroyTargetNodes();
    bam_hdr_destroy(hdr);
    sam_close(fp);
    fclose(of);

    fprintf(stderr, "Found %" PRIu32 " sights.\n", total);
    return 0;
}
