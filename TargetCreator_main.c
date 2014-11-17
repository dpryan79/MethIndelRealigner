#include "realigner.h"
//This is more of an example than anything else

void usage(char *prog) {
    fprintf(stderr, "%s file.bam output.bed\n", prog);
}

int main(int argc, char *argv[]) {
    htsFile *fp;
    bam_hdr_t *hdr;
    uint32_t total = 0;
    FILE *of;

    if(argc != 3) {
        usage(argv[0]);
        return 1;
    }

    //Open
    fp = sam_open(argv[1], "rb");
    of = fopen(argv[2], "w");
    hdr = sam_hdr_read(fp);
    initTargetNodes();
    //Iterate
    findInDels(fp, hdr, 10, 17);
    //Filter
    total = depthFilter(4);
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
