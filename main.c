#include <stdio.h>
#include <string.h>

int main_zip(int argc, char *argv[]);

static void usage(FILE *fp) {
    fprintf(fp,
"\n"
"Parfait : PARallel FAsta to 2bIt Tool in C\n"
"\n"
"Commands:\n"
"zip        FASTA file to .fz file\n"
"help       display this help message\n"
"\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) { usage(stderr); return 1; }

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { usage(stdout); return 0; }

        argv += 1;
        argc = 2;
    }

    int ret = 0;

    if (strcmp(argv[1], "zip") == 0)    ret = main_zip(argc - 1, argv + 1);
    else {
        fprintf(stderr, "unrecognized command '%s'\n", argv[1]);
        return 1;
    }

    return ret;
}