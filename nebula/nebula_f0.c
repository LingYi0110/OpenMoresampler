#include "nebula.h"

#include "../logging.h"
#include "../external/libgvps/gvps.h"
#include "../external/fastapprox/fastexp.h"
#include "../external/fastapprox/fastlog.h"
#include "../external/ciglet/ciglet.h"
#include "nebula_noise.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


#define NEB_FS    8000.0f
#define NEB_NBAND 36
#define NEB_NF0   128

// clamped 20*log10 via fast-log
static inline float neb_fast_mag2db(float x)
{
    float e = x < 1.0e-10f ? 1.0e-10f : x;
    return fastlog(e) * 8.685889988860346;
}

// F0 estimation: nebula_est, f0_from_wave, biased_f0_estimation

typedef struct {
    float sigma;
    float *ptrans;
} neb_vit_ctx;


// gvps F0-viterbi transition
static float neb_vit_ftran(void *ctx, int ds, int t)
{
    neb_vit_ctx *c = (neb_vit_ctx *)ctx;
    (void)t;
    return c->ptrans[ds];
}

static int neb_vit_fntran(void *ctx, int t)
{
    neb_vit_ctx *c = (neb_vit_ctx *)ctx;
    (void)t;
    return (int)(5.0 * c->sigma);
}

// 2-state voicing transitions
static float neb_voice_ftran(void *ctx, int ds, int t)
{
    float *p = (float *)ctx;
    (void)t;
    if (ds == 0)
        return 1.0f - *p;
    return *p;
}

static int neb_voice_fntran(void *ctx, int t)
{
    (void)ctx;
    (void)t;
    return 2;
}


// voiced/unvoiced posterior from peak log-likelihood
float *nebula_voicing_posterior(
    int Count,
    int use_uv_offset,
    const float *x,    // per-frame peak loglik
    float mu,
    float var,
    float *out_logZ)
{
    float *result = (float *)calloc((size_t)Count, sizeof(float));
    if (out_logZ)
        *out_logZ = 0.0f;
    if (Count <= 0)
        return result;

    float inv_s = 1.0 / sqrt(6.283185307179586 * var);
    float half_inv = -0.5f / var;
    float logsum = 0.0f;

    if (use_uv_offset != 0) {
        for (int i = 0; i < Count; i++) {
            float d = x[i] - mu;
            float v_pdf = fastexp(d * d * half_inv) * inv_s;
            float du = x[i] + 4.7800002f;
            float u_pdf = fastexp(du * du * -25.00000055879356f) *
                          2.8209479f;
            float den = v_pdf + u_pdf;
            result[i] = v_pdf / den;
            logsum += fastlog(den);
        }
    } else {
        for (int i = 0; i < Count; i++) {
            float d = x[i] - mu;
            float v_pdf = fastexp(d * d * half_inv) * inv_s;
            float u_pdf = fastpow2(x[i] * -147.21378f) *
                           102.04082f;
            float den = v_pdf + u_pdf;
            result[i] = v_pdf / den;
            logsum += fastlog(den);
        }
    }
    if (out_logZ)
        *out_logZ = logsum;
    return result;
}

static float neb_logp(float p)
{
    return fastlog(p);
}

static float neb_logadd(float a, float b)
{
    float hi = a > b ? a : b;
    float lo = a > b ? b : a;
    return hi + fastlog(fastexp(lo - hi) + 1.0f);
}

// 5-iter EM around the voicing posterior
static void nebula_voicing_em(const float *obs, int Count, int use_uv_offset,
                              float *mu, float *var, float *pchange)
{
    if (!obs)
        return;
    float *alpha_v = (float *)calloc((size_t)(Count > 0 ? Count : 0), sizeof(*alpha_v));
    float *alpha_u = (float *)calloc((size_t)(Count > 0 ? Count : 0), sizeof(*alpha_u));
    float *beta_v = (float *)calloc((size_t)(Count > 0 ? Count : 0), sizeof(*beta_v));
    float *beta_u = (float *)calloc((size_t)(Count > 0 ? Count : 0), sizeof(*beta_u));
    if (!alpha_v || !alpha_u || !beta_v || !beta_u) {
        free(alpha_v); free(alpha_u); free(beta_v); free(beta_u);
        return;
    }

    for (int iter = 0; iter < 5; ++iter) {
        float unused_logz = 0.0f;
        float *emission = nebula_voicing_posterior(
            Count, use_uv_offset, obs, *mu, *var, &unused_logz);
        float lstay = neb_logp(1.0f - *pchange);
        float lchange = neb_logp(*pchange);
        float sw = 0.0f;
        float swx = 0.0f;
        float swx2 = 0.0f;
        float schange = 0.0f;

        if (!emission)
            break;
        if (Count > 0) {
            alpha_v[0] = neb_logp(emission[0]);
            alpha_u[0] = neb_logp(1.0f - emission[0]);
        }
        if (Count > 1) {
            for (int t = 1; t < Count; ++t) {
                float lv = neb_logp(emission[t]);
                float lu = neb_logp(1.0f - emission[t]);
                alpha_v[t] = lv + neb_logadd(alpha_v[t - 1] + lstay,
                                             alpha_u[t - 1] + lchange);
                alpha_u[t] = lu + neb_logadd(alpha_u[t - 1] + lstay,
                                             alpha_v[t - 1] + lchange);
            }
            beta_v[Count - 1] = 0.0f;
            beta_u[Count - 1] = 0.0f;
            for (int t = Count - 2; t >= 0; --t) {
                float lv = neb_logp(emission[t + 1]);
                float lu = neb_logp(1.0f - emission[t + 1]);
                beta_v[t] = neb_logadd(lstay + lv + beta_v[t + 1],
                                        lchange + lu + beta_u[t + 1]);
                beta_u[t] = neb_logadd(lchange + lv + beta_v[t + 1],
                                        lstay + lu + beta_u[t + 1]);
            }
            float logz = neb_logadd(alpha_v[Count - 1], alpha_u[Count - 1]);
            for (int t = 1; t < Count; ++t) {
                float lv = neb_logp(emission[t]);
                float lu = neb_logp(1.0f - emission[t]);
                float g = fastexp(alpha_v[t] + beta_v[t] - logz);
                sw += g;
                swx += g * obs[t];
                swx2 += g * obs[t] * obs[t];
                schange += fastexp(alpha_v[t - 1] + lchange + lu
                                        + beta_u[t] - logz);
                schange += fastexp(alpha_u[t - 1] + lchange + lv
                                        + beta_v[t] - logz);
            }
            *pchange = schange / sw;
            *mu = swx / sw;
            *var = swx2 / sw - *mu * *mu;
        }
        free(emission);
    }

    free(alpha_v); free(alpha_u); free(beta_v); free(beta_u);
}

// per-frame likelihood energy for unvoiced observations
static float *nebula_likelihood_energy(float **mm, int Count)
{
    float *out = (float *)calloc((size_t)Count, sizeof(float));
    if (Count <= 0)
        return out;
    for (int t = 0; t < Count; t++) {
        float *row = mm[t];
        float s = 0.0f;
        for (int k = 0; k < NEB_NF0; k++) {
            float L = row[k];
            if (L > -20.0f)
                s -= L * fastexp(L);
        }
        float v = 4.8520303f - s;
        if (v < 1.0e-10f)
            v = 1.0e-10f;
        out[t] = v;
    }
    return out;
}

// free f0, aux, f_axis, lmap rows
void nebula_result_free(nebula_result *r)
{
    if (!r)
        return;
    free(r->f0);
    free(r->aux);
    free(r->f_axis);
    if (r->nfrm != 0) {
        for (int i = 0; i < r->nfrm; i++)
            free(r->lmap_rows[i]);
    }
    free(r->lmap_rows);
    free(r);
}

// resample, preprocess, IF/SNR, GMM likelihood, Viterbi, refine F0
nebula_result *nebula_est(
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
    int unused12, int unused13, int unused14,
    float *wave, int nx_in,
    float fs)
{
    (void)unused12;
    (void)unused13;
    (void)unused14;

    int nx = nx_in;
    float ratio = NEB_FS / fs;
    int nrs = nx;
    float *rs = cig_rresample(wave, nx, ratio, &nrs);
    nx = nrs;
    int Count = (int)(nx * 0.00012500001f / thop);

    float v20 = cos(fmin * 0.00025000001 * M_PI);
    float b1[2];
    b1[0] = 1.0f;
    b1[1] = -1.0f;
    float a_c[2];
    a_c[0] = 1.0f;
    a_c[1] = (float)(-((-(v20 * -2.0) -
                         sqrt(v20 * -2.0 * (v20 * -2.0) -
                              (v20 * 8.0 - 7.0) * 4.0)) * 0.5));
    float *y1 = cig_filter(b1, 2, a_c, 2, rs, nx);
    float *rev = (float *)malloc((size_t)nx * sizeof(float));
    float *y2;
    if (nx > 0) {
        for (int i = 0; i < nx; i++)
            rev[i] = y1[nx - 1 - i];
        y2 = cig_filter(b1, 2, a_c, 2, rev, nx);
        for (int i = 0; i < nx; i++)
            y1[i] = y2[nx - 1 - i];
    } else {
        y2 = cig_filter(b1, 2, a_c, 2, rev, nx);
    }
    free(rev);
    free(y2);
    free(rs);
    float *x = y1;

    float *x2 = (float *)malloc((size_t)nx * sizeof(float));
    if (nx > 0) {
        for (int i = 0; i < nx; i++)
            x2[i] = x[i] * x[i];
        float *rms = cig_moving_avg(x2, nx, 256.0f);
        for (int i = 0; i < nx; i++)
            rms[i] = sqrt(rms[i]);
        float peak = rms[0];
        for (int i = 1; i < nx; i++)
            if (rms[i] > peak)
                peak = rms[i];
        float th = peak * dither_level;
        for (int i = 0; i < nx; i++) {
            if (th > rms[i]) {
                x[i] = nebula_normal() * (th - rms[i]) + x[i];
            }
        }
        free(rms);
    } else {
        float *rms = cig_moving_avg(x2, nx, 256.0f);
        free(rms);
    }
    free(x2);

    int *naxis = (int *)calloc((size_t)Count, sizeof(int));
    for (int j = 0; j < Count; j++)
        naxis[j] = (int)(j * (thop * 8000.0));

    float *fc = (float *)malloc(NEB_NBAND * sizeof(float));
    for (int k = 0; k < NEB_NBAND; k++)
        fc[k] = k * 0.039941143f + 1.60206f;
    for (int k = 0; k < NEB_NBAND; k++)
        fc[k] = (float)(pow(10.0, fc[k]) * 0.000125);

    float *f_axis = (float *)malloc(NEB_NF0 * sizeof(float));
    for (int k = 0; k < NEB_NF0; k++)
        f_axis[k] = k * 0.011007402f + 1.60206f;
    for (int k = 0; k < NEB_NF0; k++)
        f_axis[k] = (float)pow(10.0, f_axis[k]);

    float **SNR1 = (float **)calloc((size_t)Count, sizeof(float *));
    float **SNR2 = (float **)calloc((size_t)Count, sizeof(float *));
    float **IF1 = (float **)calloc((size_t)Count, sizeof(float *));
    float **IF2 = (float **)calloc((size_t)Count, sizeof(float *));
    if (Count > 0) {
        for (int t = 0; t < Count; t++) {
            SNR1[t] = (float *)calloc(NEB_NBAND, sizeof(float));
            SNR2[t] = (float *)calloc(NEB_NBAND, sizeof(float));
            IF1[t] = (float *)calloc(NEB_NBAND, sizeof(float));
            IF2[t] = (float *)calloc(NEB_NBAND, sizeof(float));
        }
    }
    estimate_ifsnr(x, nx, naxis, Count, fc, fc, NEB_NBAND,
                   "nuttall98", "hanning",
                   SNR1, IF1, SNR2, IF2);

    for (int t = 0; t < Count; t++) {
        for (int b = 0; b < NEB_NBAND; b++) {
            SNR1[t][b] = neb_fast_mag2db(SNR1[t][b]);
            SNR2[t][b] = neb_fast_mag2db(SNR2[t][b]);
        }
    }

    neb_gmm *bands[NEB_NBAND];
    for (int b = 0; b < NEB_NBAND; b++)
        bands[b] = neb_gmm_load_band(b);

    float **mm = (float **)calloc((size_t)Count, sizeof(*mm));
    for (int t = 0; t < Count; ++t)
        mm[t] = (float *)calloc(NEB_NF0, sizeof(*mm[t]));
    neb_make_lmap(bands, SNR1, SNR2, IF1, IF2, Count, f_axis, mm);

    for (int b = 0; b < NEB_NBAND; b++)
        neb_gmm_free(bands[b]);

    if (Count != 0) {
        for (int t = 0; t < Count; t++)
            free(SNR1[t]);
        free(SNR1);
        for (int t = 0; t < Count; t++)
            free(SNR2[t]);
        free(SNR2);
        for (int t = 0; t < Count; t++)
            free(IF1[t]);
        free(IF1);
        for (int t = 0; t < Count; t++)
            free(IF2[t]);
        free(IF2);
    } else {
        free(SNR1);
        free(SNR2);
        free(IF1);
        free(IF2);
    }

    if (lmap_smooth > 0.0f) {
        float *col = (float *)calloc((size_t)Count, sizeof(float));
        for (int k = 0; k < NEB_NF0; k++) {
            float halford = lmap_smooth / f_axis[k] * (2.0f / thop);
            int w = (int)halford;
            for (int t = 0; t < Count; t++)
                col[t] = fastexp(mm[t][k]);
            float *sm = cig_moving_avg(col, Count, halford * 0.5f);
            if (Count - w > w) {
                for (int t = w; t < Count - w; t++)
                    mm[t][k] = fastlog(sm[t]);
            }
            free(sm);
        }
        free(col);
    }

    int i1 = 0;
    for (i1 = 0; i1 != NEB_NF0; ++i1)
        if (fmin <= f_axis[i1])
            break;
    int i2 = 0;
    for (i2 = 0; i2 != NEB_NF0; ++i2)
        if (fmax <= f_axis[i2])
            break;
    int nstate = i2 - i1;

    float blockb = viterbi_sigma / (log10(fmax) - log10(fmin)) * (thop * 128.0f);
    float *ptrans = (float *)calloc(NEB_NF0, sizeof(float));
    {
        double gauss_n = 0.3989422804014327 / blockb;
        float sum = 0.0;
        double exp2v = 1.0;
        for (int i = 0; i < NEB_NF0; i++) {
            float v = exp2v * gauss_n;
            ptrans[i] = v;
            sum += v;
            if (i != NEB_NF0 - 1) {
                double t = i + 1;
                exp2v = exp(t * t * (-0.5 / (blockb * blockb)));
            }
        }
        float inv_sum = 1.0 / sum;
        for (int i = 0; i < NEB_NF0; i++)
            ptrans[i] *= inv_sum;
    }

    gvps_obsrv *obs = gvps_obsrv_create_full(Count, nstate);
    if (Count > 0 && i1 < i2) {
        for (int t = 0; t < Count; t++) {
            for (int i = 0; i < nstate; i++) {
                float L = mm[t][i1 + i];
                obs->slice[t]->pair[i].state = i;
                obs->slice[t]->pair[i].p = fastexp(L);
            }
        }
    }
    neb_vit_ctx vctx;
    vctx.sigma = blockb;
    vctx.ptrans = ptrans;
    int *path = (int *)malloc((size_t)Count * sizeof(int));
    gvps_full_sampled_static(path, &vctx, obs, neb_vit_ftran, neb_vit_fntran);
    if (Count > 0) {
        for (int t = 0; t < Count; t++)
            path[t] += i1;
    }
    gvps_obsrv_free(obs);
    free(ptrans);

    float *f0 = (float *)calloc((size_t)Count, sizeof(float));
    for (int t = 0; t < Count; t++) {
        int pk = path[t];
        float idx = (float)pk;
        if (pk != 127 && pk > 0)
            cig_qifft(mm[t], pk, &idx);
        int i0 = (int)floor(idx);
        if (i0 == 127)
            f0[t] = f_axis[127];
        else
            f0[t] = (idx - i0) * (f_axis[i0 + 1] - f_axis[i0])
                    + f_axis[i0];
    }

    if (skip_voicing == 0) {
        const int peak_lik = use_peak_likelihood != 0;
        float *voice_obs;
        float pchange = thop / voicing_time_const;
        float mu = peak_lik ? -2.0f : 1.0f;
        float var = 1.0f;
        float final_logz = 0.0f;

        if (peak_lik) {
            voice_obs = (float *)calloc((size_t)Count, sizeof(*voice_obs));
            for (int t = 0; t < Count; ++t) {
                float v = mm[t][path[t]];
                if (v < -4.7800002f)
                    v = -4.7800002f;
                voice_obs[t] = v;
            }
        } else {
            voice_obs = nebula_likelihood_energy(mm, Count);
        }

        nebula_voicing_em(voice_obs, Count, peak_lik,
                          &mu, &var, &pchange);
        float final_mu = voicing_mu == 0.0f ? mu : voicing_mu;
        float final_var = voicing_var == 0.0f ? var : voicing_var;
        float *post = nebula_voicing_posterior(Count, peak_lik,
                                         voice_obs, final_mu, final_var,
                                         &final_logz);
        gvps_obsrv *voice_vit = gvps_obsrv_create_full(Count, 2);
        int *voice_path = (int *)malloc((size_t)Count * sizeof(*voice_path));
        if (Count > 0 && post && voice_vit) {
            for (int t = 0; t < Count; ++t) {
                voice_vit->slice[t]->pair[0].state = 0;
                voice_vit->slice[t]->pair[0].p = 1.0f - post[t];
                voice_vit->slice[t]->pair[1].state = 1;
                voice_vit->slice[t]->pair[1].p = post[t];
            }
        }
        gvps_full_sampled_static(voice_path, &pchange, voice_vit,
                          neb_voice_ftran, neb_voice_fntran);
        if (Count > 0) {
            for (int t = 0; t < Count; ++t)
                f0[t] *= (float)voice_path[t];
        }
        free(voice_path);
        if (voice_vit)
            gvps_obsrv_free(voice_vit);
        free(post);
        free(voice_obs);
    } else if (Count > 0) {
        int *voice_path = (int *)malloc((size_t)Count * sizeof(*voice_path));
        for (int t = 0; t < Count; ++t)
            voice_path[t] = 1;
        for (int t = 0; t < Count; ++t)
            f0[t] *= (float)voice_path[t];
        free(voice_path);
    }

    free(path);
    free(naxis);
    free(fc);
    free(x);

    nebula_result *r = (nebula_result *)malloc(sizeof(*r));
    r->f0 = f0;
    r->lmap_rows = mm;
    r->aux = nebula_likelihood_energy(mm, Count);
    r->f_axis = f_axis;
    r->nfrm = Count;
    r->lmap_width = NEB_NF0;
    return r;
}

// recursive 6-point search on kalman L; returns best log10 bias
static float refine_range_search(
    float *z, float *Rmeas, float *Rbase, int Count,
    float lo, float hi, int depth)
{
    if (depth == 0)
        return (hi + lo) * 0.5f;

    float *Qtry = (float *)calloc((size_t)Count, sizeof(float));
    float *Ptmp = (float *)calloc((size_t)Count, sizeof(float));
    float step = (hi - lo) * 0.16666667f; // 1/6
    float Lbest = -1.0e10f;
    int best_i = 1;

    for (int i = 1; i != 6; i++) {
        float scale = pow(10.0f, (float)i * step + lo);
        float Lcur = 0.0f;
        for (int t = 0; t < Count; t++)
            Qtry[t] = Rbase[t] * scale;
        float *xf = cig_kalmanf1d(z, Qtry, Rmeas, Count, 1.0f, Ptmp, &Lcur);
        free(xf);
        if (Lcur > Lbest) {
            Lbest = Lcur;
            best_i = i;
        }
    }
    free(Qtry);
    free(Ptmp);

    float lo2 = step * (float)(best_i - 1) + lo;
    float hi2 = (float)(best_i + 1) * step + lo;
    return refine_range_search(z, Rmeas, Rbase, Count, lo2, hi2, depth - 1);
}

// interpolate voiced F0 control points, then moving-average of 5
static void gapfill_smooth(float *f0, int nf0)
{
    int nc = 1;
    float *cx = (float *)calloc((size_t)nf0 + 3, sizeof(float));
    float *cy = (float *)calloc((size_t)nf0 + 3, sizeof(float));
    float *unit = (float *)malloc((size_t)nf0 * sizeof(float));
    for (int k = 0; k < nf0; k++) {
        unit[k] = (float)k;
        if (f0[k] != 0.0f) {
            cx[nc] = (float)k;
            cy[nc] = f0[k];
            nc++;
        }
    }
    if (nc == 1) {
        free(cx); free(cy); free(unit);
        return;
    }
    cx[0] = 0.0f;
    cy[0] = cy[1];
    cx[nc] = (float)nf0;
    cy[nc] = cy[nc - 1];
    float *itp = cig_interp(cx, cy, nc + 1, unit, nf0);
    float *sm = cig_moving_avg(itp, nf0, 5.0f);
    for (int k = 0; k < nf0; k++)
        f0[k] = sm[k];
    free(cx); free(cy); free(unit); free(itp); free(sm);
}

// optional Kalman / harmonic refinement of an F0 contour
float *biased_f0_estimation(
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
    float *wave, int nx,
    float fs,
    float *f0_in, int nf0)
{
    float hop = thop;

    (void)fmin;
    (void)fmax;
    (void)use_peak_likelihood;
    (void)skip_voicing;
    (void)dither_level;
    (void)voicing_time_const;
    (void)viterbi_sigma;
    (void)lmap_smooth;
    (void)voicing_mu;
    (void)voicing_var;

    if (nf0 <= 0)
        return (float *)calloc(1, sizeof(float));

    {
        int all_uv = 1;
        for (int k = 0; k < nf0; k++) {
            if (f0_in[k] > 0.0f) {
                all_uv = 0;
                break;
            }
        }
        if (all_uv)
            return (float *)calloc((size_t)nf0, sizeof(float));
    }

    int *naxis = (int *)calloc((size_t)nx, sizeof(int));
    float *f0 = (float *)calloc((size_t)nf0, sizeof(float));
    for (int k = 0; k < nf0; k++) {
        float t = (float)k * (hop * fs);
        naxis[k] = (int)roundf(t);
        f0[k] = f0_in[k] * (1.0f / fs); // ratio track
    }

    gapfill_smooth(f0, nf0);

    for (int pass = 0; pass < n_pass; pass++) {
        int decf = (int)roundf(64.0f / (float)nhar);
        float *blk;
        float *phase;
        int ncyc = 0;

        if (decf < 1)
            decf = 1;

        // gap-fill again into blk, then *2
        {
            float *cx = (float *)calloc((size_t)nf0 + 3, sizeof(float));
            float *cy = (float *)calloc((size_t)nf0 + 3, sizeof(float));
            float *unit = (float *)malloc((size_t)nf0 * sizeof(float));
            int nc = 1;
            for (int k = 0; k < nf0; k++) {
                unit[k] = (float)k;
                if (f0[k] != 0.0f) {
                    cx[nc] = (float)k;
                    cy[nc] = f0[k];
                    nc++;
                }
            }
            if (nc == 1) {
                free(cx);
                free(cy);
                free(unit);
                break;
            }
            cx[0] = 0.0f;
            cy[0] = cy[1];
            cx[nc] = (float)nf0;
            cy[nc] = cy[nc - 1];
            blk = cig_interp(cx, cy, nc + 1, unit, nf0);
            free(cx);
            free(cy);
            free(unit);
        }
        for (int k = 0; k < nf0; k++)
            blk[k] += blk[k];

        // phase integration along the F0 contour
        {
            float cap = (float)naxis[nf0 - 1];
            int maxph = naxis[nf0 - 1] > 0 ? naxis[nf0 - 1] : 1;
            float ph = 0.0f;
            float p = 0.5f;
            float blk_cur = blk[0] > 1e-12f ? blk[0] : 1e-12f;
            int cur = 0;
            phase = (float *)calloc((size_t)maxph + 8, sizeof(float));
            ncyc = 0;
            if (cap > 0.0f) {
                for (;;) {
                    float cand = p / blk_cur + ph;
                    float interp, next;
                    if (cand > (float)naxis[cur]) {
                        while (++cur != nf0) {
                            if (cand <= (float)naxis[cur])
                                break;
                        }
                        if (cur >= nf0)
                            cur = nf0 - 1;
                        else
                            cur -= 1;
                    } else if (cur != 0) {
                        cur -= 1;
                    }
                    if (cur >= nf0 - 1) {
                        interp = blk[nf0 - 1] > 1e-12f ? blk[nf0 - 1] : 1e-12f;
                        next = ph + 1.0f / interp;
                        if (ncyc < maxph + 4)
                            phase[ncyc++] = next;
                        if (cap <= next)
                            break;
                        ph = next;
                        blk_cur = interp;
                        continue;
                    }
                    {
                        float denom =
                            (float)(naxis[cur + 1] - naxis[cur]);
                        float frac =
                            denom != 0.0f
                                ? (cand - (float)naxis[cur]) / denom
                                : 0.0f;
                        interp = (1.0f - frac) * blk[cur]
                                 + frac * blk[cur + 1];
                        if (interp < 1e-12f)
                            interp = 1e-12f;
                        next = ph + 1.0f / interp;
                        if (ncyc < maxph + 4)
                            phase[ncyc++] = next;
                        if (cap <= next)
                            break;
                        ph = next;
                        blk_cur = interp;
                    }
                    if (ncyc > maxph + 2)
                        break;
                }
            }
        }

        float *cyc_axis = (float *)calloc(ncyc > 0 ? (size_t)ncyc : 1u, sizeof(float));
        for (int k = 0; k < ncyc; k++)
            cyc_axis[k] = roundf((float)(k << 8) * 0.5f);
        int v210 = (ncyc > 0) ? (int)cyc_axis[ncyc - 1] : 0;
        if (v210 < 1)
            v210 = 1;

        float *res_axis = (float *)malloc((size_t)v210 * sizeof(float));
        for (int k = 0; k < v210; k++)
            res_axis[k] = (float)k;
        float *wave_axis = (float *)malloc((size_t)(nx > 0 ? nx : 1) * sizeof(float));
        for (int k = 0; k < nx; k++)
            wave_axis[k] = (float)k;
        float *naxisf = (float *)calloc((size_t)nf0, sizeof(float));
        for (int k = 0; k < nf0; k++)
            naxisf[k] = (float)naxis[k];

        float *phase2samp = cig_interp(cyc_axis, phase, ncyc > 0 ? ncyc : 1,
                                res_axis, v210);
        float *samp2phase = cig_interp(phase, cyc_axis, ncyc > 0 ? ncyc : 1,
                                naxisf, nf0);
        float *res_wave = cig_interp(wave_axis, wave, nx, phase2samp, v210);

        free(naxisf);
        free(wave_axis);
        free(res_axis);
        free(phase);
        free(cyc_axis);
        free(blk);

        float dec_bw = (float)decf * 0.00390625f; // decf/256
        float cutoff = (1.0f / (float)decf) * 0.5f;
        float *fir = cig_winfir(32, cutoff, 0.0f, (char *)"lowpass", (char *)"hamming");
        float *lp = cig_convolution(res_wave, fir, v210, 32);

        int v95 = v210 / decf;
        if (v95 < 1)
            v95 = 1;
        float *dec = (float *)calloc((size_t)v95, sizeof(float));
        for (int j = 0; j < v95; j++) {
            int idx = 16 + j * decf; // conv length v210+32-1, still indexed from +16
            dec[j] = lp[idx];
        }
        free(lp);
        free(fir);

        int *wf0 = (int *)calloc((size_t)nf0, sizeof(int));
        for (int k = 0; k < nf0; k++)
            wf0[k] = (int)roundf(samp2phase[k] * (1.0f / (float)decf));

        float *hidx = (float *)malloc((size_t)nhar * sizeof(float));
        float *bw = (float *)malloc((size_t)nhar * sizeof(float));
        for (int h = 0; h < nhar; h++) {
            float top = (float)nhar * dec_bw;
            hidx[h] = (float)h * ((top - dec_bw) / (float)(nhar - 1))
                      + dec_bw;
            bw[h] = dec_bw;
        }

        float **ha = (float **)calloc((size_t)nf0, sizeof(float *));
        float **period_out = (float **)calloc((size_t)nf0, sizeof(float *));
        for (int f = 0; f < nf0; f++) {
            ha[f] = (float *)calloc((size_t)nhar, sizeof(float));
            period_out[f] = (float *)calloc((size_t)nhar, sizeof(float));
        }

        // SNR into ha, IF into period_out on the decimated wave (nuttall83)
        estimate_ifsnr(dec, v95, wf0, nf0, hidx, bw, nhar,
                       "nuttall83", "nuttall83",
                       ha, period_out, NULL, NULL);

        float **hf = (float **)calloc((size_t)nf0, sizeof(float *));
        for (int f = 0; f < nf0; f++)
            hf[f] = (float *)calloc((size_t)nhar, sizeof(float));

        for (int f = 0; f < nf0; f++) {
            for (int h = 0; h < nhar; h++) {
                double mag = ha[f][h];
                if (mag < 1.0e-10)
                    mag = 1.0e-10;
                double hnr_db = 20.0 * log10(mag);
                float conf;
                if (hnr_db > -70.0) {
                    float inner = fastexp(hnr_db * 0.1286 + 6.6856);
                    conf = fastexp(-1.755 - 1.0 / (inner + 0.38207));
                } else {
                    conf = fastexp(hnr_db * 0.058142 + 0.20593);
                }
                int untrusted = (conf > 0.1)
                            || (fabsf(period_out[f][h] - hidx[h])
                                > 0.3 * bw[h]);
                if (untrusted) {
                    float t = conf * 100.0f;
                    hf[f][h] = t * t;
                } else {
                    float t = conf / (float)(h + 1);
                    hf[f][h] = t * t;
                }
            }
        }
        for (int f = 0; f < nf0; f++)
            for (int h = 0; h < nhar; h++)
                period_out[f][h] /= ((float)h + 1.0f) * dec_bw;

        // inverse-variance fuse; conf copied into ha, weights into hf
        float *fused = (float *)calloc((size_t)nf0, sizeof(float));
        float *fvar = (float *)calloc((size_t)nf0, sizeof(float));
        for (int f = 0; f < nf0; f++) {
            float wsum = 0.0f;
            float acc = 0.0f, vacc = 0.0f;
            for (int h = 0; h < nhar; h++)
                ha[f][h] = hf[f][h];
            for (int h = 0; h < nhar; h++) {
                float conf = ha[f][h];
                float w = 0.5f / conf;
                hf[f][h] = w;
                wsum += w;
            }
            if (wsum > 0.0f) {
                for (int h = 0; h < nhar; h++)
                    hf[f][h] /= wsum;
            }
            for (int h = 0; h < nhar; h++) {
                float w = hf[f][h];
                float conf = ha[f][h];
                acc += period_out[f][h] * w;
                vacc += w * w * conf; // weight^2 * original hf conf
            }
            fused[f] = acc;
            fvar[f] = vacc;
        }

        for (int f = 0; f < nf0; f++) {
            free(ha[f]);
            free(period_out[f]);
            free(hf[f]);
        }
        free(ha);
        free(period_out);
        free(hf);
        free(hidx);
        free(bw);
        free(wf0);
        free(dec);
        free(phase2samp);
        free(samp2phase);
        free(res_wave);

        float *R;
        {
            float half = (float)kalman_win; // may be 0
            float *dev = cig_moving_avg(fused, nf0, half);
            for (int k = 0; k < nf0; k++) {
                float d = fused[k] - dev[k];
                dev[k] = d * d;
            }
            R = cig_moving_avg(dev, nf0, half);
            free(dev);
        }
        float log_bias = refine_range_search(fused, fvar, R, nf0, -3.0f, 0.0f, 5);
        float kbias = pow(10.0f, log_bias);
        float *Q = (float *)calloc((size_t)nf0, sizeof(float));
        for (int k = 0; k < nf0; k++)
            Q[k] = R[k] * kbias;
        free(R);

        float *P = (float *)calloc((size_t)nf0, sizeof(float));
        float *kf = cig_kalmanf1d(fused, Q, fvar, nf0, fused[0], P, NULL);
        float *ks = cig_kalmans1d(kf, P, Q, nf0);
        for (int k = 0; k < nf0; k++)
            f0[k] = ks[k] * f0[k];

        free(fused);
        free(fvar);
        free(Q);
        free(P);
        free(kf);
        free(ks);
    }

    for (int k = 0; k < nf0; k++) {
        if (f0_in[k] != 0.0f)
            f0[k] = f0[k] * fs;
        else
            f0[k] = 0.0f;
    }

    free(naxis);
    return f0;
}
