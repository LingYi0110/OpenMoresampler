#include "nebula.h"
#include "../logging.h"
#include "../external/ciglet/ciglet.h"
#include "nebula_noise.h"

#include <math.h>
#include <stdlib.h>

// two-pass F0: coarse nebula_est, then refine
float *f0_from_wave(
    float fmin, float fmax,
    float thop,
    int use_peak_likelihood,
    int skip_voicing,
    float dither_level,
    float voicing_time_const,
    float viterbi_sigma,
    float lmap_smooth,
    float voicing_mu,
    float voicing_var,
    int n_pass, int nhar, int kalman_win,
    float *x, int nx,
    int fs,
    int *out_nfrm)
{
    srand(0);

    float fs_f = (float)fs;

    nebula_result *r0 = nebula_est(
        fmin, fmax,
        0.03f,
        use_peak_likelihood, skip_voicing,
        dither_level, voicing_time_const, viterbi_sigma,
        lmap_smooth, voicing_mu, voicing_var,
        n_pass, nhar, kalman_win,
        x, nx, fs_f);

    int nfrm = r0->nfrm;
    int width = r0->lmap_width;

    float *row_max = (float *)calloc((size_t)nfrm, sizeof(float));
    for (int i = 0; i < nfrm; i++)
        row_max[i] = maxfp(r0->lmap_rows[i], width);

    int *idx = (int *)malloc((size_t)nfrm * sizeof(int));
    float *sorted = cig_sort(row_max, nfrm, idx);

    int ntop = nfrm <= 10 ? nfrm : 10;
    // overwrite argsort values with the selected f_axis peaks
    for (int t = 0; t < ntop; t++) {
        int fi = idx[nfrm - 1 - t];
        int best_k = 0;
        float best_v = -100.0f;
        if (width > 0) {
            for (int k = 0; k < width; k++) {
                float v = r0->lmap_rows[fi][k];
                if (v > best_v) {
                    best_v = v;
                    best_k = k;
                }
            }
        }
        sorted[t] = r0->f_axis[best_k];
    }

    float median = cig_medianfp(sorted, ntop);
    log_tee_printf("estimate_f0_range: median = %f Hz\n", (double)median);

    free(row_max);
    free(idx);
    // argsort value buffer is left allocated
    nebula_result_free(r0);

    float *xd = (float *)calloc((size_t)nx, sizeof(float));
    if (nx > 0) {
        for (int i = 0; i < nx; i++)
            xd[i] = fabsf(x[i]);
        float peak = maxfp(xd, nx);

        for (int i = 0; i < nx; i++) {
            float xi = x[i];
            xd[i] = nebula_normal() * peak * 0.01 + xi;
        }
    }

    nebula_result *r1 = nebula_est(
        fmin, fmax,
        thop,
        1, skip_voicing,
        dither_level, voicing_time_const, viterbi_sigma,
        2.0f, voicing_mu, voicing_var,
        n_pass, nhar, kalman_win,
        xd, nx, fs_f);

    *out_nfrm = r1->nfrm;
    free(xd);

    int all_uv = 1;
    if (r1->nfrm > 0) {
        for (int i = 0; i < r1->nfrm; i++) {
            if (r1->f0[i] != 0.0f) {
                all_uv = 0;
                break;
            }
        }
    } else {
        all_uv = 1;
    }

    if (all_uv || r1->nfrm <= 0) {
        log_tee_printf(
            "biased_f0_estimation: input is unvoiced; no need for refinement.\n");
        float *out = (float *)calloc((size_t)r1->nfrm, sizeof(float));
        nebula_result_free(r1);
        return out;
    }

    log_tee_printf("biased_f0_estimation: refining...\n");
    {
        float *out = biased_f0_estimation(
            fmin, fmax, thop,
            1, skip_voicing,
            dither_level, voicing_time_const, viterbi_sigma,
            2.0f, voicing_mu, voicing_var,
            n_pass, nhar, kalman_win,
            x, nx, fs_f,
            r1->f0, r1->nfrm);
        nebula_result_free(r1);
        log_tee_printf("biased_f0_estimation: completed.\n");
        return out;
    }
}
