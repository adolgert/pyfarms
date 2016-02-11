#include <stdio.h>
#include <stdlib.h>
#include "herd_size_factor.h"


int main(char** argv, int argc) {
    int pr_idx;
    double* size_factor=0;
    int herds[]={1, 3, 13, 7, 10, 7, 5, 2, 10, 12, 2510, 100530, 285310, 567290, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000, 36000};
    int herd_cnt=sizeof(herds)/sizeof(int);
    size_factor=malloc(herd_cnt*sizeof(double));

    herd_size_factor(herds, herd_cnt, size_factor);

    for (pr_idx=0; pr_idx<herd_cnt; ++pr_idx) {
        printf("%g, ", pr_idx+1, size_factor[pr_idx]);
    }

    free(size_factor);

    return 0;
}