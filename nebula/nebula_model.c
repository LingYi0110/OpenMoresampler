#include "nebula.h"

#include "../external/fastapprox/fastexp.h"
#include "../external/fastapprox/fastlog.h"
#include "../external/ciglet/ciglet.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


#define NEB_NBAND     36
#define NEB_NF0       128
#define NEB_NMIX      16
#define NEB_NORD      6
#define NEB_NFEAT     5
#define NEB_GMM_BYTES 2752

#include "model/nebula_lcal.inc"
#include "model/nebula_gmm.inc"


// free one band GMM
void neb_gmm_free(neb_gmm *g)
{
    if (!g)
        return;
    if (g->mu) {
        for (int i = 0; i < g->nmix; i++)
            free(g->mu[i]);
        free(g->mu);
    }
    if (g->sigma) {
        for (int i = 0; i < g->nmix; i++)
            free(g->sigma[i]);
        free(g->sigma);
    }
    if (g->inv11) {
        for (int i = 0; i < g->nmix; i++)
            free(g->inv11[i]);
        free(g->inv11);
    }
    if (g->s21s11i) {
        for (int i = 0; i < g->nmix; i++)
            free(g->s21s11i[i]);
        free(g->s21s11i);
    }
    free(g->w);
    free(g->det11);
    free(g->condsigma);
    free(g);
}

// load one band GMM from the embedded blob; LU of the 5x5 feature block
neb_gmm *neb_gmm_load_band(int band)
{
    if (band < 0 || band >= NEB_NBAND)
        return NULL;
    const float *src = (const float *)(NEBULA_GMM_BLOB + band * NEB_GMM_BYTES);
    neb_gmm *g = (neb_gmm *)calloc(1, sizeof(*g));

    g->nmix = NEB_NMIX;
    g->nord = NEB_NORD;
    g->w = (float *)malloc(NEB_NMIX * sizeof(float));
    g->mu = (float **)calloc(NEB_NMIX, sizeof(float *));
    g->sigma = (float **)calloc(NEB_NMIX, sizeof(float *));
    g->det11 = (float *)malloc(NEB_NMIX * sizeof(float));
    g->inv11 = (float **)calloc(NEB_NMIX, sizeof(float *));
    g->s21s11i = (float **)calloc(NEB_NMIX, sizeof(float *));
    g->condsigma = (float *)malloc(NEB_NMIX * sizeof(float));

    memcpy(g->w, src, NEB_NMIX * sizeof(float));
    const float *p = src + NEB_NMIX;

    for (int m = 0; m < NEB_NMIX; m++) {
        g->mu[m] = (float *)malloc(NEB_NORD * sizeof(float));
        g->sigma[m] = (float *)malloc(NEB_NORD * NEB_NORD * sizeof(float));
        memcpy(g->mu[m], p, NEB_NORD * sizeof(float));
        p += NEB_NORD;
        memcpy(g->sigma[m], p, NEB_NORD * NEB_NORD * sizeof(float));
        p += NEB_NORD * NEB_NORD;
    }

    for (int m = 0; m < NEB_NMIX; m++) {
        float *sig = g->sigma[m];
        float *A = (float *)calloc(NEB_NFEAT * NEB_NFEAT, sizeof(float));

        for (int r = 0; r < NEB_NFEAT; r++)
            for (int c = 0; c < NEB_NFEAT; c++)
                A[r + c * NEB_NFEAT] = sig[r * NEB_NORD + c];
        cig_lu(A, NEB_NFEAT);
        g->inv11[m] = A;

        float det = A[0];
        for (int k = 1; k < NEB_NFEAT; k++)
            det *= A[k * NEB_NFEAT + k];
        g->det11[m] = logf(det);

        float *rhs = (float *)malloc(NEB_NFEAT * sizeof(float));
        for (int k = 0; k < NEB_NFEAT; k++)
            rhs[k] = sig[k * NEB_NORD + 5];
        cig_lusolve(A, rhs, NEB_NFEAT);
        g->s21s11i[m] = rhs;

        float cs = sig[35];
        for (int k = 0; k < NEB_NFEAT; k++)
            cs -= sig[30 + k] * rhs[k];
        g->condsigma[m] = cs;
    }
    g->ready = 1;
    return g;
}

// Conditional GMM log density for one band.
static float *neb_condgmm_density(const neb_gmm *g, const float *feat,
                                  const float *axis, int width)
{
    int nmix = g->nmix;
    int nfeat = g->nord - 1;
    float *weights = (float *)calloc((size_t)nmix, sizeof(float));
    float *means = (float *)calloc((size_t)nmix, sizeof(float));
    float *sigma = (float *)calloc((size_t)nmix, sizeof(float));
    float *logw = (float *)calloc((size_t)nmix, sizeof(float));
    float *y = (float *)calloc((size_t)(2 * g->nord), sizeof(float));
    double normalizer = nfeat * 1.837877066409345;
    float maximum = -1.0e10f;

    for (int i = 0; i < nmix; ++i) {
        const float *mu = g->mu[i];
        float mean = mu[nfeat];
        means[i] = mean;
        sigma[i] = g->condsigma[i] < 4.0 ? 4.0 : g->condsigma[i];
        for (int k = 0; k < nfeat; ++k) {
            y[k] = feat[k] - mu[k];
            mean += y[k] * g->s21s11i[i][k];
        }
        means[i] = mean;
        cig_lusolve(g->inv11[i], y, nfeat);
        float quad = 0.0f;
        for (int k = 0; k < nfeat; ++k)
            quad += (feat[k] - mu[k]) * y[k];
        float value = fastlog(g->w[i]) - normalizer
                    + (quad + g->det11[i]) * -0.5;
        logw[i] = value;
        if (maximum < value)
            maximum = value;
    }
    float sum = 0.0f;
    for (int i = 0; i < nmix; ++i) {
        float value = logw[i] - maximum;
        if (value < -10.0)
            value = -10.0;
        weights[i] = fastexp(value);
        sum += weights[i];
    }
    float invsum = 1.0 / sum;
    for (int i = 0; i < nmix; ++i)
        weights[i] *= invsum;

    float *density = (float *)calloc((size_t)width, sizeof(float));
    for (int i = 0; i < nmix; ++i) {
        if (weights[i] < 1.0e-5)
            continue;
        float scale = weights[i] / sqrt(6.283185307179586 * sigma[i]);
        float invsigma = 1.0 / sigma[i];
        for (int k = 0; k < width; ++k) {
            float d = axis[k] - means[i];
            float q = d * d * invsigma;
            float value = q <= 50.0 ? fasterexp(-0.5 * q) : 0.0f;
            density[k] += value * scale;
        }
    }
    for (int k = 0; k < width; ++k)
        density[k] = fastlog(density[k]);
    free(weights);
    free(means);
    free(sigma);
    free(logw);
    free(y);
    return density;
}

static void neb_make_lmap_impl(
    neb_gmm **bands, float **SNR1, float **SNR2, float **IF1, float **IF2,
    int nfrm, int nband, const float *axis, int width, float **mm)
{
#pragma omp parallel for
    for (int t = 0; t < nfrm; ++t) {
        float *row = mm[t];
        memset(row, 0, (size_t)width * sizeof(float));
        for (int b = 0; b < nband; ++b) {
            float feat[5];
            feat[0] = SNR1[t][b];
            feat[1] = SNR1[t][b <= 6 ? b + 7 : b - 7];
            feat[2] = SNR2[t][b];
            feat[3] = IF1[t][b] * 8000.0;
            feat[4] = 8000.0 * IF2[t][b];
            float *density = neb_condgmm_density(bands[b], feat, axis, width);
            for (int k = 0; k < width; ++k)
                row[k] += density[k];
            free(density);
        }
        float invb = 1.0 / nband;
        float maximum = -1.0e10f;
        for (int k = 0; k < width; ++k) {
            float value = row[k] * invb - NEBULA_LCAL_MEAN[k];
            row[k] = value;
            if (maximum < value)
                maximum = value;
        }
        float sum = 0.0f;
        for (int k = 0; k < width; ++k) {
            row[k] -= maximum;
            sum += fastexp(row[k]);
        }
        float logsum = fastlog(sum);
        for (int k = 0; k < width; ++k)
            row[k] -= logsum;
    }
}

void neb_make_lmap(neb_gmm **bands, float **SNR1, float **SNR2,
                   float **IF1, float **IF2, int nfrm,
                   const float *f_axis, float **mm)
{
    neb_make_lmap_impl(bands, SNR1, SNR2, IF1, IF2,
                       nfrm, NEB_NBAND, f_axis, NEB_NF0, mm);
}
