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

// accumulate one band of the conditional GMM onto out_ll
static void neb_condgmm_accum(
    const neb_gmm *g,
    const float feat[5],
    const float *f_axis,
    float *out_ll)
{
    float f0_mu[NEB_NMIX];
    float f0_sigma[NEB_NMIX];
    float f0_lL[NEB_NMIX];
    float f0_w[NEB_NMIX];
    float dx[NEB_NFEAT];
    float density[NEB_NF0];
    float lwsum = -1.0e10f;
    const float v103 = NEB_NFEAT * 1.837877066409345f;
    float wsum = 0.0f;

    memset(density, 0, sizeof(density));

    for (int i = 0; i < g->nmix; i++) {
        float *mu = g->mu[i];
        float mu0 = mu[NEB_NFEAT];
        float cs = g->condsigma[i];
        float y[NEB_NFEAT];
        float quad = 0.0f;

        if (cs < 4.0f)
            cs = 4.0f;
        f0_sigma[i] = cs;

        for (int k = 0; k < NEB_NFEAT; k++) {
            dx[k] = feat[k] - mu[k];
            mu0 += g->s21s11i[i][k] * dx[k];
        }
        f0_mu[i] = mu0;

        memcpy(y, dx, sizeof(y));
        cig_lusolve(g->inv11[i], y, NEB_NFEAT);
        for (int k = 0; k < NEB_NFEAT; k++)
            quad += dx[k] * y[k];

        float lw = fastlog(g->w[i]) - v103
           + (quad + g->det11[i]) * -0.5f;
        f0_lL[i] = lw;
        if (lw > lwsum)
            lwsum = lw;
    }

    for (int i = 0; i < g->nmix; i++) {
        float lw = f0_lL[i] - lwsum;
        if (lw < -10.0f)
            lw = -10.0f;
        f0_w[i] = fastexp(lw);
        wsum += f0_w[i];
    }
    for (int i = 0; i < g->nmix; i++)
        f0_w[i] /= wsum;

    for (int i = 0; i < g->nmix; i++) {
        if (f0_w[i] < 1.0e-5f)
            continue;
        float invs = 1.0f / f0_sigma[i];
        float scale = f0_w[i] / sqrtf(6.283185307179586f * f0_sigma[i]);
        for (int k = 0; k < NEB_NF0; k++) {
            float d = f_axis[k] - f0_mu[i];
            float qq = d * d * invs;
            if (qq > 50.0f)
                continue;
            density[k] += fasterexp(-0.5f * qq) * scale;
        }
    }
    for (int k = 0; k < NEB_NF0; ++k)
        out_ll[k] += fastlog(density[k]);
}

// per-frame GMM likelihood map, subtract Lcal, log-softmax
void neb_make_lmap(
    neb_gmm **bands,
    float **SNR1, float **SNR2, float **IF1, float **IF2,
    int nfrm,
    const float *f_axis,
    float **mm)
{
#pragma omp parallel for
    for (int t = 0; t < nfrm; t++) {
        float acc[NEB_NF0];
        float invb = 1.0f / NEB_NBAND;
        float vmax = -1.0e10f;
        float z = 0.0f;

        memset(acc, 0, sizeof(acc));
        for (int b = 0; b < NEB_NBAND; b++) {
            float feat[5];
            feat[0] = SNR1[t][b];
            feat[1] = SNR1[t][b <= 6 ? b + 7 : b - 7];
            feat[2] = SNR2[t][b];
            feat[3] = IF1[t][b] * 8000.0f;
            feat[4] = 8000.0f * IF2[t][b];
            neb_condgmm_accum(bands[b], feat, f_axis, acc);
        }
        for (int k = 0; k < NEB_NF0; k++) {
            float v = invb * acc[k] - NEBULA_LCAL_MEAN[k];
            mm[t][k] = v;
            if (v > vmax)
                vmax = v;
        }
        for (int k = 0; k < NEB_NF0; k++) {
            mm[t][k] -= vmax;
            z += fastexp(mm[t][k]);
        }
        float logz = fastlog(z);
        for (int k = 0; k < NEB_NF0; k++)
            mm[t][k] -= logz;
    }
}
