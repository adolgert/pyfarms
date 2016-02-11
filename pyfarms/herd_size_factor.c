/*! The herd size factor is described generally in documentation.
    This is an attempt to replicate the NAADSM code in order to figure
    out what the factors really are.
 */
#include <stdio.h>
#include "gsl/gsl_histogram.h"


void herd_size_factor(int* herds, int herd_cnt, double* size_factor) {
    gsl_histogram* histogram=0;
    gsl_histogram_pdf* histogram_pdf=0;
    int max_size=-1;
    int max_idx, incr_idx, sf_idx;
    double first_x, last_x;
    double area;
    size_t bin;
    double lower, upper;
    double frac;

    // Following airborne-spread-model.c
    for (max_idx=0; max_idx<herd_cnt; ++max_idx) {
        max_size=(herds[max_idx]>max_size) ? herds[max_idx] : max_size;
    }

    histogram=gsl_histogram_alloc(max_size);
    gsl_histogram_set_ranges_uniform(histogram, 0.5, (double) max_size+0.5);
    for (incr_idx=0; incr_idx<herd_cnt; ++incr_idx) {
        gsl_histogram_increment(histogram, (double) herds[incr_idx]);
    }

    // Here following prob_dist.c to make the pdf.
    gsl_histogram_scale(histogram, 1.0/gsl_histogram_sum(histogram));
    histogram_pdf=gsl_histogram_pdf_alloc(gsl_histogram_bins(histogram));
    gsl_histogram_pdf_init(histogram_pdf, histogram);
    first_x=gsl_histogram_min(histogram);
    last_x=gsl_histogram_max(histogram);

    for (sf_idx=0; sf_idx<herd_cnt; ++sf_idx) {
        if (herds[sf_idx]<first_x) {
            area=0.0;
        } else if (herds[sf_idx]>=last_x) {
            area=1.0;
        } else {
            gsl_histogram_find(histogram, herds[sf_idx], &bin);
            gsl_histogram_get_range(histogram, bin, &lower, &upper);
        }
        frac=(herds[sf_idx]-lower)/(upper-lower);
        area=histogram_pdf->sum[bin]+frac*(histogram_pdf->sum[bin+1]-histogram_pdf->sum[bin]);
        size_factor[sf_idx]=2*area;
    }
    gsl_histogram_pdf_free(histogram_pdf);
    gsl_histogram_free(histogram);
}
