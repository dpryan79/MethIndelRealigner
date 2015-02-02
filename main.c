#include <stdio.h>
#include <string.h>
#include "htslib/hts.h"
#include "version.h" //This has the VERSION define

int realigner_main(int argc, char *argv[]);
int TargetCreator_main(int argc, char *argv[]);

void usage_main() {
    fprintf(stderr, "MethIndelRealigner: A tool for local realignment of bisulfite sequencing data around InDels.\n");
    fprintf(stderr, "Version: %s (using HTSlib version %s)\n", VERSION, hts_version());
    fprintf(stderr, "Usage: MethIndelRealigner <command> [options] <file(s)>\n");
    fprintf(stderr, "\n"
"Commands:\n"
"    TargetCreator  Determine target regions for local realignment\n"
"    realign        Perform local realignment around target InDels.\n"
);
}

int main(int argc, char *argv[]) {
    if(argc == 1) {
        usage_main();
        return 0;
    } else if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        usage_main();
        return 0;
    } else if(strcmp(argv[1], "TargetCreator") == 0) {
        return TargetCreator_main(argc-1, argv+1);
    } else if(strcmp(argv[1], "realign") == 0) {
        return realigner_main(argc-1, argv+1);
    } else {
        fprintf(stderr, "Unknown command!\n");
        usage_main();
        return -1;
    }
}

