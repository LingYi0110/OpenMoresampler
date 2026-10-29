#include "lmix.h"
#include "lmix_utils.h"

#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#include "../config.h"
#include "../external/ciglet/ciglet.h"
#include "lmix_lf.h"
#include "../external/libllsm/envelope.h"
#include "../logging.h"
#include "../resampler.h"

// 1 - exp((Te-T0)*x) - Ta*x
static FP_TYPE lmix_lf_efunc(FP_TYPE x, void *env)
{
    lmix_lfparam *p = (lmix_lfparam *)env;
    return 1.0 - exp_3((p->Te - p->T0) * x) - p->Ta * x;
}

// (x^2+wg^2)*sin*A + sin*x + wg*exp(-x*Te) - wg*cos
static FP_TYPE lmix_lf_afunc_cb(FP_TYPE x, void *env)
{
    lmix_lfparam *p = (lmix_lfparam *)env;
    FP_TYPE weighted_sine = p->sin_wgTe * p->A;
    return (x * x + p->wg * p->wg) * weighted_sine +
           p->sin_wgTe * x +
           p->wg * exp_3(-x * p->Te) -
           p->wg * p->cos_wgTe;
}

// Brent root finder for the LF model equations.
static FP_TYPE lmix_lf_brent(FP_TYPE (*func)(FP_TYPE, void *), FP_TYPE xmin, FP_TYPE xmax, void* env) {
  const FP_TYPE eps = 1e-8;
  FP_TYPE a, b, c, s, d = 0.0f;
  FP_TYPE fa, fb, fc, fs;
  a = xmin; b = xmax;
  fa = func(a, env); fb = func(b, env);
  if(fa * fb >= 0) return (a + b) / 2.0;
  if(fabs(fa) < fabs(fb)) {
    c = b; b = a; a = c;
    fc = fb; fb = fa; fa = fc;
  } else {
    c = a;
    fc = fa;
  }
  int mflag = 1;
  int niter = 0;
  while(fabs(fb) > eps && fabs(a - b) > eps) {
    if(fa != fc && fb != fc) {
      s = a * fb * fc / (fa - fb) / (fa - fc)
        + b * fa * fc / (fb - fa) / (fb - fc)
        + c * fa * fb / (fc - fa) / (fc - fb); // Inv. Quad Interpolation
    } else {
      s = b - fb * (b - a) / (fb - fa); // Secant
    }
    if((s < (3.0 * a + b) / 4.0 || s > b) ||
       (mflag   && fabs(s - b) >= fabs(b - c) * 0.5) ||
       (! mflag && fabs(s - b) >= fabs(c - d) * 0.5) ||
       (mflag   && fabs(b - c) < eps) ||
       (! mflag && fabs(c - d) < eps)) {
      s = (a + b) / 2.0;
      if(s == b || s == a) break; // maximum precision reached
      mflag = 1;
    } else {
      mflag = 0;
    }
    fs = func(s, env);
    d = c; c = b;
    fc = fb;
    if(fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if(fabs(fa) < fabs(fb)) {
      FP_TYPE x = b; b = a; a = x;
      x = fb; fb = fa; fa = x;
    }
    niter ++;
  }
  return b;
}

lmix_lfparam lmix_lfparam_from_model(lfmodel model)
{
    lmix_lfparam ret;
    ret.T0 = model.T0;
    ret.Te = model.T0 * model.te;
    ret.Tp = model.T0 * model.tp;
    ret.Ta = model.T0 * model.ta;
    ret.a = 0.0f;
    ret.wg = M_PI / ret.Tp;
    ret.sin_wgTe = sin_3(ret.wg * ret.Te);
    ret.cos_wgTe = cos_3(ret.wg * ret.Te);

    FP_TYPE e = lmix_lf_brent(lmix_lf_efunc, 1.0f, 2.0f / (ret.Ta + 1e-9), &ret);
    FP_TYPE e_Te_T0 = exp_3(e * (ret.Te - ret.T0));
    ret.A = (1.0 - e_Te_T0) / (e * e * ret.Ta) +
            (ret.Te - ret.T0) * e_Te_T0 / (e * ret.Ta);
    ret.e = e;

    ret.a = lmix_lf_brent(lmix_lf_afunc_cb, 0.0f, 1.0e9f, &ret);

    ret.E0 = -model.Ee / (exp_3(ret.a * ret.Te) * ret.sin_wgTe);
    return ret;
}

// LF * lip radiation * last-band vt_magn; optional minphase into phse_out
static FP_TYPE *lmix_frame_get_harmonic_amp(
    lmix_frame *frame, FP_TYPE *phse_out)
{
    lfmodel lf = lfmodel_from_rd(frame->rd, 1.0f / frame->f0, 1.0f);
    int nhar = frame->nhar;
    FP_TYPE *freq = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
    if (nhar > 0) {
        for (int i = 0; i != nhar; ++i)
            freq[i] = ((FP_TYPE)i + 1.0f) * frame->f0;
    }
    FP_TYPE *lf_phse = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
    FP_TYPE *lf_amp = lmix_lfmodel_spectrum(lf, freq, nhar, lf_phse);
    FP_TYPE *lip_phse = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
    FP_TYPE *lip_amp = llsm_liprad(freq, nhar, 1.5f, lip_phse);
    FP_TYPE *vt_lin;

    if (nhar <= 1) {
        if (nhar > 0)
            lf_amp[0] = 1.0f;
        vt_lin = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
        if (nhar <= 0)
            goto finish;
    } else {
        for (int i = 1; i != nhar; ++i)
            lf_amp[i] = lf_amp[i] / (((FP_TYPE)i + 1.0f) * lf_amp[0]);
        lf_amp[0] = 1.0f;
        vt_lin = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
    }

    {
        int stride = frame->vt_stride;
        FP_TYPE *src = frame->vt_magn + stride - 1;
        for (int i = 0; i < nhar; ++i) {
            FP_TYPE v = *src;
            if (v < -100.0f)
                v = -100.0f;
            vt_lin[i] = expf(v);
            src += stride;
        }
    }

finish:
    {
        FP_TYPE *vt_phse = llsm_harmonic_minphase(vt_lin, nhar);
        if (nhar > 0) {
            if (phse_out != NULL) {
                for (int i = 0; i != nhar; ++i) {
                    vt_lin[i] = vt_lin[i] * lf_amp[i] * lip_amp[i];
                    phse_out[i] = vt_phse[i] + lip_phse[i] +
                        frame->vs_har_phse[i];
                }
            } else {
                for (int i = 0; i != nhar; ++i)
                    vt_lin[i] = vt_lin[i] * lf_amp[i] * lip_amp[i];
            }
        }
        free(vt_phse);
    }
    free(freq);
    free(lip_phse);
    free(lip_amp);
    free(lf_phse);
    free(lf_amp);
    return vt_lin;
}

llsm_layer0 *lmix_to_layer0(lmix_data *lmix)
{
    llsm_layer0 *layer0 = (llsm_layer0 *)calloc(1, sizeof(llsm_layer0));
    layer0->conf.nfrm = lmix->conf.nfrm;
    layer0->conf.nhop = lmix->conf.nhop;
    layer0->conf.nhar = lmix->conf.nhar;
    layer0->conf.nhare = lmix->conf.nhare;
    layer0->conf.nnos = lmix->conf.nnos;
    layer0->conf.noswarp = lmix->conf.noswarp;
    layer0->conf.mvf = lmix->conf.mvf;
    layer0->conf.nosf = lmix->conf.nosf;
    layer0->conf.thop = lmix->conf.thop;
    layer0->conf.nnosband = lmix->conf.nnosband;

    int nband = lmix->conf.nnosband - 1;
    layer0->conf.nosbandf = (FP_TYPE *)calloc((size_t)nband, sizeof(FP_TYPE));
    if (nband > 0) {
        int i = 0;
        do {
            layer0->conf.nosbandf[i] = lmix->conf.nosbandf[i];
            ++i;
        } while (i != nband);
    }

    if (lmix->conf.nfrm <= 0) {
        layer0->frames = NULL;
        return layer0;
    }

    layer0->frames = (llsm_frame **)calloc(
        (size_t)lmix->conf.nfrm, sizeof(llsm_frame *));

    for (int i = 0; i < lmix->conf.nfrm; ++i) {
        lmix_frame *src = lmix->frames[i];
        llsm_nosframe *noise = src->noise;
        llsm_frame *dst = llsm_create_frame(
            src->nhar, lmix->conf.nhare, noise->nnos, noise->nchannel);
        layer0->frames[i] = dst;
        llsm_copy_nosframe(dst->noise, src->noise);

        if (src->wave != NULL) {
            if (dst->wave == NULL) {
                dst->wave = llsm_create_wave(
                    src->wave->samples,
                    src->wave->nsample,
                    src->wave->fs);
            }
            llsm_copy_wave(dst->wave, src->wave);
        }

        dst->f0 = src->f0;
        if (src->f0 > 0.0f) {
            FP_TYPE *amp = lmix_frame_get_harmonic_amp(
                src, dst->sinu->phse);
            int nhar = dst->sinu->nhar;
            if (nhar > 0) {
                for (int k = 0; k != nhar; ++k)
                    dst->sinu->ampl[k] = amp[k];
            }
            free(amp);

            if (dst->noise->nchannel > 0) {
                for (int ch = 0; ch < dst->noise->nchannel; ++ch) {
                    llsm_sinframe *eenv = dst->noise->eenv[ch];
                    int enh = eenv->nhar;
                    FP_TYPE ph0 = dst->sinu->phse[0];
                    if (enh <= 0)
                        continue;
                    for (int k = 0; k != enh; ++k)
                        eenv->phse[k] += ((FP_TYPE)k + 1.0f) * ph0;
                }
            }
        }
    }
    return layer0;
}

// 16-byte blob from shift_spectrum; caller frees delta then the blob
typedef struct lmix_shift_blob {
    int nnos;
    FP_TYPE energy; // 0.5 * ln(E_new / E_old), or copied from anchor
    FP_TYPE loud;   // ln(loudness gain), or copied from anchor
    FP_TYPE *delta; // residual geometric envelope
} lmix_shift_blob;

// uniform Rd axis, 76 points, 0.02 .. 3.02 step 0.04
static const FP_TYPE lmix_rd_axis[76] = {
    0.02f, 0.06f, 0.10f, 0.14f, 0.18f, 0.22f, 0.26f, 0.30f,
    0.34f, 0.38f, 0.42f, 0.46f, 0.50f, 0.54f, 0.58f, 0.62f,
    0.66f, 0.70f, 0.74f, 0.78f, 0.82f, 0.86f, 0.90f, 0.94f,
    0.98f, 1.02f, 1.06f, 1.10f, 1.14f, 1.18f, 1.22f, 1.26f,
    1.30f, 1.34f, 1.38f, 1.42f, 1.46f, 1.50f, 1.54f, 1.58f,
    
    1.62f, 1.66f, 1.70f, 1.74f, 1.78f, 1.82f, 1.86f, 1.90f,
    1.94f, 1.98f, 2.02f, 2.06f, 2.10f, 2.14f, 2.18f, 2.22f,
    2.26f, 2.30f, 2.34f, 2.38f, 2.42f, 2.46f, 2.50f, 2.54f,
    2.58f, 2.62f, 2.66f, 2.70f, 2.74f, 2.78f, 2.82f, 2.86f,
    2.90f, 2.94f, 2.98f, 3.02f
};

// nonlinear mapped Rd, 76 points
static const FP_TYPE lmix_rd_mapped[76] = {
    0.02f, 0.17746f, 0.43753f, 0.66692f, 0.84729f, 0.98347f,
    1.08087f, 1.17392f, 1.25649f, 1.32908f, 1.40173f, 1.47488f,
    1.54063f, 1.59777f, 1.65342f, 1.71218f, 1.77129f, 1.82458f,
    1.86939f, 1.90908f, 1.94913f, 1.99283f, 2.03990f, 2.08734f,
    2.13110f, 2.16819f, 2.19832f, 2.22371f, 2.24745f, 2.27216f,
    2.29938f, 2.32964f, 2.36271f, 2.39785f, 2.43399f, 2.46993f,
    2.50445f, 2.53652f, 2.56546f, 2.59096f, 2.61312f, 2.63236f,
    2.64923f, 2.66434f, 2.67822f, 2.69131f, 2.70397f, 2.71641f,
    2.72881f, 2.74126f, 2.75381f, 2.76648f, 2.77930f, 2.79226f,
    2.80536f, 2.81859f, 2.83193f, 2.84541f, 2.85898f, 2.87267f,
    2.88645f, 2.90032f, 2.91425f, 2.92820f, 2.94213f, 2.95597f,
    2.96960f, 2.98181f, 2.98587f, 2.99078f, 2.99628f, 3.00205f,
    3.00774f, 3.01291f, 3.01715f, 3.02f
};

static void lmix_shift_blob_free(void *p)
{
    lmix_shift_blob *b = (lmix_shift_blob *)p;

    if (b == NULL)
        return;
    free(b->delta);
    free(b);
}

// ITU-R 468-ish weighting used by both shift functions
static FP_TYPE lmix_loudness_weight(FP_TYPE freq)
{
    FP_TYPE f = freq < 300.0f ? 300.0f : freq;
    FP_TYPE f2 = f * f;
    FP_TYPE h = f2 * f2 * (f2 + 1440000.0f)
        / ((f2 + 160000.0f) * (f2 + 160000.0f) * (f2 + 9610000.0f));
    return (FP_TYPE)pow((double)h, 0.33);
}

// map rd through the nonlinear table, subtract 0.01*Mt, invert, clamp to [0.02, 3.02]
static FP_TYPE lmix_shift_rd(FP_TYPE rd, FP_TYPE mt)
{
    FP_TYPE xi1 = lmix_rd_axis[75] + lmix_rd_axis[75] - lmix_rd_axis[74];
    FP_TYPE *tmp = cig_interpu(lmix_rd_axis[0], xi1, (FP_TYPE *)lmix_rd_mapped,
                      76, &rd, 1);
    FP_TYPE mapped = tmp[0] - 0.01 * mt;
    free(tmp);
    tmp = cig_interp((FP_TYPE *)lmix_rd_mapped, (FP_TYPE *)lmix_rd_axis,
                     76, &mapped, 1);
    FP_TYPE out = tmp[0];
    free(tmp);
    if (out > 3.02f)
        return 3.02f;
    if (out < 0.02f)
        return 0.02f;
    return out;
}

// windowed-sinc upsample of vt_magn onto the FFT axis, then floor near/above harmonic Nyquist
static FP_TYPE *lmix_vt_to_faxis(
    FP_TYPE *vt_freq, FP_TYPE *vt_magn, int nvt,
    FP_TYPE *faxis, int nbin, FP_TYPE fhar)
{
    int n4 = 4 * nvt + 4;
    FP_TYPE *freq4 = (FP_TYPE *)calloc((size_t)n4, sizeof(FP_TYPE));
    FP_TYPE *conv = (FP_TYPE *)calloc((size_t)n4, sizeof(FP_TYPE));
    if (nvt > 0) {
        freq4[0] = 0.0f;
        freq4[1] = 0.25f * vt_freq[0];
        freq4[2] = 0.50f * vt_freq[0];
        freq4[3] = 0.75f * vt_freq[0];
        for (int i = 1; i < nvt; ++i) {
            FP_TYPE prev = vt_freq[i - 1];
            FP_TYPE df = vt_freq[i] - prev;
            freq4[4 * i] = prev;
            freq4[4 * i + 1] = prev + 0.25f * df;
            freq4[4 * i + 2] = prev + 0.50f * df;
            freq4[4 * i + 3] = prev + 0.75f * df;
        }
    }

    FP_TYPE *kernel = (FP_TYPE *)calloc(41, sizeof(FP_TYPE));
    for (int i = 0; i < 41; ++i)
        kernel[i] = 0.54 - 0.46 * cos_2(i * M_PI / 20.0);
    for (int i = -20; i < 0; ++i) {
        FP_TYPE x = i * M_PI / 4.0;
        kernel[i + 20] *= sin_2(x) / x;
    }
    for (int i = 1; i < 20; ++i) {
        FP_TYPE x = i * M_PI / 4.0;
        kernel[i + 20] *= sin_2(x) / x;
    }

    int k = n4 <= 20 ? n4 : 20;
    if (k > 0 && nvt > 0) {
        for (int i = 0; i < k; ++i)
            conv[i] = kernel[i + 20] * vt_magn[0];
    }
    if (nvt > 0) {
        int pos = 0;
        for (int h = 1; h <= nvt; ++h) {
            int end = nvt > h + 3 ? 4 * h + 20 : n4;
            while (pos < end) {
                conv[pos] += kernel[pos + 20 - 4 * h] * vt_magn[h - 1];
                ++pos;
            }
            if (h == nvt)
                break;
            pos = h <= 3 ? 0 : 4 * h - 16;
        }
    }
    free(kernel);

    FP_TYPE *windowed = cig_interp(freq4, conv, n4, faxis, nbin);
    free(conv);
    free(freq4);
    FP_TYPE fnyq_ext = vt_freq[nvt - 1] + vt_freq[nvt - 1] - vt_freq[nvt - 2];
    FP_TYPE *direct = cig_interpu(vt_freq[0], fnyq_ext, vt_magn, nvt, faxis, nbin);

    FP_TYPE fade0 = 0.690779983997345f;
    for (int i = 0; i < nbin; ++i) {
        FP_TYPE fade;
        if (faxis[i] >= fhar)
            fade = 0.0f;
        else if (faxis[i] >= fhar - 1000.0f)
            fade = (fhar - faxis[i]) * 0.000690779983997345;
        else
            fade = fade0;
        FP_TYPE floorv = direct[i] + fade;
        if (windowed[i] > floorv)
            windowed[i] = floorv;
    }
    free(direct);
    return windowed;
}

// gender formant warp (g) plus ME/Mo/Mr spectral filters
static FP_TYPE *lmix_gender_warp(
    FP_TYPE *spec, int nbin, FP_TYPE *faxis,
    FP_TYPE fhar, FP_TYPE f0_new, resampler_args *args)
{
    FP_TYPE fnyq = faxis[nbin - 1];
    FP_TYPE src_x[5];
    FP_TYPE dst_x[5];
    src_x[0] = 0.0f;
    src_x[1] = f0_new;
    src_x[2] = 0.8f * fhar;
    if (src_x[2] > 4000.0f)
        src_x[2] = 4000.0f;
    src_x[3] = 0.9f * fhar;
    if (src_x[3] > 8000.0f)
        src_x[3] = 8000.0f;
    src_x[4] = fnyq;
    dst_x[0] = 0.0f;
    dst_x[1] = f0_new;
    dst_x[2] = src_x[2] * (1.0f - 0.01f * args->g);
    if (src_x[3] - 1000.0f < dst_x[2])
        dst_x[2] = src_x[3] - 1000.0f;
    if (dst_x[2] < 100.0f)
        dst_x[2] = 100.0f;
    dst_x[3] = src_x[3];
    dst_x[4] = src_x[4];

    FP_TYPE *map = cig_interp(dst_x, src_x, 5, faxis, nbin);
    FP_TYPE *out = cig_interpu(faxis[0], fnyq + fnyq - faxis[nbin - 2],
                      spec, nbin, map, nbin);
    free(map);

    if (args->ME != 0.0f && nbin > 0) {
        FP_TYPE me = args->ME;
        FP_TYPE k = me >= 0.0f ? 0.005f : 0.0025f;
        for (int i = 0; i < nbin; ++i)
            out[i] = out[i] + me * k * out[i];
    }
    if (args->Mo != 0.0f && nbin > 0) {
        FP_TYPE two_nyq = fnyq + fnyq;
        FP_TYPE a = exp_2(-600.0f * (FP_TYPE)M_PI / two_nyq);
        FP_TYPE b = -2.0f * exp_2(-300.0f * (FP_TYPE)M_PI / two_nyq);
        FP_TYPE c = 1.0f - (b - a);
        FP_TYPE href;
        FP_TYPE omega0 = f0_new / fnyq * (FP_TYPE)(2.0 * M_PI);

        {
            FP_TYPE e = 1.0f - ((-a) * 0.99998659f + (-0.99998641f) * b);
            FP_TYPE d = a * 1.8702625e-7f;
            FP_TYPE denom = e * e + d * d;

            href = (FP_TYPE)sqrt((double)((e * c / denom) * (e * c / denom)
                + (d * c / denom) * (d * c / denom)));
        }
        for (int i = 0; i < nbin; ++i) {
            FP_TYPE theta = (0.5f / fnyq) * faxis[i] * (FP_TYPE)(2.0 * M_PI)
                - omega0 + (FP_TYPE)M_PI;
            FP_TYPE s1 = sin_2(theta);
            FP_TYPE c1 = cos_2(theta);
            FP_TYPE num_re = 1.0f - c1 * b - (-a) * s1;
            FP_TYPE num_im = sin_2(theta + theta) * (-a) + c1 * b;
            FP_TYPE denom = num_re * num_re + num_im * num_im;
            FP_TYPE mag = (FP_TYPE)sqrt((double)((num_re * c / denom) * (num_re * c / denom)
                + (num_im * c / denom) * (num_im * c / denom))) / href;
            out[i] -= log_2(mag) * args->Mo * 0.003333333333333334f;
        }
    }
    if (args->Mr != 0.0f && nbin > 0) {
        FP_TYPE mr = args->Mr;
        FP_TYPE inv = 1.0f / (fnyq + fnyq);
        FP_TYPE a = exp_2(-6283.185307179586f * inv);
        FP_TYPE b = -2.0f * exp_2(-3141.592653589793f * inv);
        FP_TYPE c = 1.0f - ((-a) + b);
        FP_TYPE omega_off = inv * 20106.19298297468f;
        FP_TYPE e = 1.0f - ((-a) * 0.99998659f + (-0.99998641f) * b);
        FP_TYPE d = a * 1.8702625e-7f;
        FP_TYPE denom0 = e * e + d * d;
        FP_TYPE href = (FP_TYPE)sqrt((double)((e * c / denom0) * (e * c / denom0)
            + (c * d / denom0) * (c * d / denom0)));
        for (int i = 0; i < nbin; ++i) {
            FP_TYPE theta = (0.5f / fnyq) * faxis[i] * (FP_TYPE)(2.0 * M_PI)
                - omega_off + (FP_TYPE)M_PI;
            FP_TYPE s1 = sin_2(theta);
            FP_TYPE c1 = cos_2(theta);
            FP_TYPE num_re = 1.0f - c1 * b - (-a) * s1;
            FP_TYPE num_im = sin_2(theta + theta) * (-a) + c1 * b;
            FP_TYPE denom = num_re * num_re + num_im * num_im;
            FP_TYPE mag = (FP_TYPE)sqrt((double)((num_re * c / denom) * (num_re * c / denom)
                + (c * num_im / denom) * (c * num_im / denom))) / href;
            FP_TYPE mr_i = mr;
            if (mr_i <= 0.0f)
                mr_i = mr_i * 0.09900990099009901f;
            FP_TYPE gain = mr_i * (mag * 0.1f) + 1.0f;
            out[i] += log_2(gain);
        }
    }
    return out;
}

// shared body of shift_spectrum (anchor == NULL) and shift_spectrum_edge
static lmix_shift_blob *lmix_shift_spectrum_core(
    lmix_frame *frame, FP_TYPE *freqwarp, resampler_args *args,
    FP_TYPE f0_new, lmix_shift_blob *anchor)
{
    FP_TYPE f0_old = frame->f0;
    if (f0_old <= 0.0f)
        return NULL;

    lmix_shift_blob *blob = (lmix_shift_blob *)calloc(1, sizeof(*blob));
    int nhar_old = frame->nhar;
    int stride = frame->vt_stride;
    FP_TYPE fhar = (FP_TYPE)nhar_old * f0_old;
    int nhar_new = (int)ceil((double)fhar / (double)f0_new);
    int nhar_max = nhar_old < nhar_new ? nhar_new : nhar_old;
    int nvt = nhar_max * stride;
    int nnos = frame->noise->nnos;
    int nfft = g_lmix_nfft;
    int nbin = nfft / 2 + 1;
    int fs_fft = 2 * (int)(freqwarp[nnos - 1] + 1000.0f);
    blob->nnos = nnos;

    FP_TYPE *faxis = llsm_uniform_faxis(nfft, fs_fft);
    FP_TYPE *vt_freq = (FP_TYPE *)calloc((size_t)nvt, sizeof(FP_TYPE));
    FP_TYPE *har_freq = (FP_TYPE *)calloc((size_t)nhar_max, sizeof(FP_TYPE));
    FP_TYPE *amp = (FP_TYPE *)calloc((size_t)nhar_max, sizeof(FP_TYPE));
    if (nvt > 0) {
        for (int i = 0; i != nvt; ++i)
            vt_freq[i] = ((FP_TYPE)i + 1.0f) * frame->vt_freq_step;
    }
    if (nhar_max > 0) {
        for (int i = 0; i != nhar_max; ++i)
            har_freq[i] = ((FP_TYPE)i + 1.0f) * f0_old;
    }
    FP_TYPE fnyq = faxis[nbin - 1] + faxis[nbin - 1] - faxis[nbin - 2];

    FP_TYPE rd_new = lmix_shift_rd(frame->rd, args->Mt);
    lfmodel lf_s = lfmodel_from_rd(frame->rd, 1.0f / f0_old, 1.0f);
    lfmodel lf_d = lfmodel_from_rd(rd_new, 1.0f / f0_old, 1.0f);
    FP_TYPE *lf_phse = (FP_TYPE *)calloc((size_t)nhar_max, sizeof(FP_TYPE));
    FP_TYPE *lf_src = lmix_lfmodel_spectrum(lf_s, har_freq, nhar_max, lf_phse);
    FP_TYPE *lf_dst = lmix_lfmodel_spectrum(lf_d, har_freq, nhar_max, NULL);
    FP_TYPE *lip = llsm_liprad(har_freq, nhar_old, 1.5f, NULL);
    FP_TYPE *spec = lmix_vt_to_faxis(
        vt_freq, frame->vt_magn, nhar_old * stride, faxis, nbin, fhar);

    if (nhar_max > 1) {
        for (int i = 1; i != nhar_max; ++i) {
            FP_TYPE d = ((FP_TYPE)i + 1.0f);
            lf_src[i] = lf_src[i] / (d * lf_src[0]);
            lf_dst[i] = lf_dst[i] / (d * lf_dst[0]);
        }
    }
    lf_src[0] = 1.0f;
    lf_dst[0] = 1.0f;
    frame->rd = rd_new;

    FP_TYPE *lf_src_f = cig_interpu(
        har_freq[0],
        har_freq[nhar_max - 1] + har_freq[nhar_max - 1] - har_freq[nhar_max - 2],
        lf_src, nhar_max, faxis, nbin);
    FP_TYPE *env_old = llsm_geometric_envelope(spec, nfft, fs_fft, freqwarp, nnos);

    FP_TYPE loud_old = 0.0f;
    FP_TYPE energy_old = 0.0f;
    if (nhar_old > 0) {
        FP_TYPE *src = frame->vt_magn + stride - 1;
        for (int i = 0; i < nhar_old; ++i) {
            amp[i] = lf_src[i] * lip[i] * exp_2(src[0]);
            src += stride;
        }
        for (int i = 0; i < nhar_old; ++i) {
            FP_TYPE w = lmix_loudness_weight(((FP_TYPE)i + 1.0f) * f0_old);
            loud_old += w * (amp[i] * amp[i]);
            energy_old += amp[i] * amp[i];
        }
    }
    free(lip);

    int skip = 0;
    FP_TYPE *warped;
    if (f0_new > f0_old) {
        FP_TYPE fcut = f0_new + 500.0f;
        FP_TYPE sx[4];
        FP_TYPE dx[4];
        sx[0] = 0.0f;
        sx[1] = f0_old;
        sx[2] = fcut;
        sx[3] = fnyq;
        dx[0] = 0.0f;
        dx[1] = f0_new;
        dx[2] = fcut;
        dx[3] = fnyq;
        FP_TYPE *map = cig_interp(dx, sx, 4, faxis, nbin);
        FP_TYPE *stretched = cig_interpu(faxis[0], fnyq, spec, nbin, map, nbin);
        if (fcut > freqwarp[0]) {
            do
                ++skip;
            while (fcut > freqwarp[skip]);
        }
        free(map);
        free(spec);
        warped = lmix_gender_warp(stretched, nbin, faxis, fhar, f0_new, args);
        free(stretched);
    } else {
        warped = lmix_gender_warp(spec, nbin, faxis, fhar, f0_new, args);
        free(spec);
    }

    FP_TYPE *env_new = llsm_geometric_envelope(warped, nfft, fs_fft, freqwarp, nnos);
    blob->delta = env_new;
    if (nnos > skip) {
        FP_TYPE *nspec = frame->noise->spec;

        if (anchor != NULL) {
            for (int i = skip; i < nnos; ++i) {
                FP_TYPE d = anchor->delta[i];

                env_new[i] = d;
                if (d > 1.0f)
                    d = 1.0f;
                nspec[i] += d;
            }
        } else {
            for (int i = skip; i < nnos; ++i) {
                FP_TYPE d = env_new[i] - env_old[i];

                env_new[i] = d;
                if (d > 1.0f)
                    d = 1.0f;
                nspec[i] += d;
            }
        }
    }
    if (skip != 0)
        memset(env_new, 0, (size_t)skip * sizeof(FP_TYPE));
    free(env_old);

    if (nhar_new > nhar_old) {
        frame->vs_har_phse = (FP_TYPE *)realloc(
            frame->vs_har_phse, (size_t)nhar_new * sizeof(FP_TYPE));
        for (int i = nhar_old; i < nhar_new; ++i)
            frame->vs_har_phse[i] = lf_phse[i];
    }

    frame->f0 = f0_new;
    frame->nhar = nhar_new;
    frame->vt_freq_step = f0_new / (FP_TYPE)stride;
    if (nvt > 0) {
        for (int i = 0; i != nhar_new * stride; ++i)
            vt_freq[i] = ((FP_TYPE)i + 1.0f) * frame->vt_freq_step;
    }
    if (nhar_new > 0) {
        for (int i = 0; i != nhar_new; ++i)
            har_freq[i] = ((FP_TYPE)i + 1.0f) * f0_new;
    }

    {
        FP_TYPE dens = log(f0_old / f0_new);

        for (int i = 0; i < nbin; ++i)
            warped[i] += dens;
    }

    FP_TYPE *lf_dst_f = cig_interpu(
        har_freq[0],
        har_freq[nhar_max - 1] + har_freq[nhar_max - 1] - har_freq[nhar_max - 2],
        lf_dst, nhar_max, faxis, nbin);
    for (int i = 0; i < nbin; ++i) {
        FP_TYPE fade;
        if (faxis[i] < fhar - 2500.0f)
            fade = 0.0f;
        else if (fhar <= faxis[i])
            fade = 1.0f;
        else
            fade = (faxis[i] - fhar + 2500.0f) * 0.0004f;
        FP_TYPE ratio = lf_src_f[i] / lf_dst_f[i];
        FP_TYPE lratio = log_2(ratio);
        FP_TYPE mixed = (1.0f - 0.01f * args->Mm) * lratio;
        FP_TYPE clamped = fade * lratio;
        if (clamped > 0.8f)
            clamped = 0.8f;
        if (clamped < -0.8f)
            clamped = -0.8f;
        if (clamped > 0.0f && mixed > clamped)
            warped[i] += mixed;
        else if (clamped <= 0.0f && clamped > mixed)
            warped[i] += mixed;
        else
            warped[i] += clamped;
    }
    free(lf_src_f);
    free(lf_dst_f);

    FP_TYPE *lip2 = llsm_liprad(har_freq, nhar_new, 1.5f, NULL);
    FP_TYPE loud_new = 0.0f;
    FP_TYPE energy_new = 0.0f;
    {
        FP_TYPE *spec_at_har = cig_interpu(
            faxis[0], fnyq, warped, nbin, har_freq, nhar_new);
        if (nhar_new > 0) {
            for (int i = 0; i < nhar_new; ++i)
                amp[i] = lf_dst[i] * lip2[i] * exp_2(spec_at_har[i]);
            for (int i = 0; i < nhar_new; ++i) {
                FP_TYPE w = lmix_loudness_weight(((FP_TYPE)i + 1.0f) * f0_new);
                loud_new += w * (amp[i] * amp[i]);
                energy_new += amp[i] * amp[i];
            }
        }
        free(spec_at_har);
        free(lip2);
    }

    if (g_synthesis_loudness_preservation != 0) {
        FP_TYPE gain;
        if (loud_new <= 0.0f)
            gain = 2.0f;
        else
            gain = (FP_TYPE)sqrt((double)(loud_old / loud_new));
        if (gain > 2.0f)
            gain = 2.0f;
        if (gain < 0.5f)
            gain = 0.5f;
        FP_TYPE gain_log = (anchor != NULL) ? anchor->loud : log(gain);
        blob->loud = gain_log;
        for (int i = 0; i < nbin; ++i)
            warped[i] += gain_log;
        if (nhar_new > 0) {
            for (int i = 0; i < nhar_new; ++i)
                amp[i] *= gain;
        }
    }

    free(frame->vt_magn);
    frame->vt_magn = cig_interpu(
        faxis[0], fnyq, warped, nbin, vt_freq, nhar_new * stride);

    energy_new = 0.0f;
    if (nhar_new > 0) {
        for (int i = 0; i < nhar_new; ++i)
            energy_new += amp[i] * amp[i];
    }
    if (anchor != NULL)
        blob->energy = anchor->energy;
    else
        blob->energy = log(energy_new / energy_old) * 0.5f;

    {
        FP_TYPE mb = args->Mb;
        FP_TYPE scale;
        if (mb < 0.0f)
            scale = (mb + 100.0f) * 0.008f + 0.2000000001f;
        else
            scale = 0.01f * mb + 1.0000000001f;
        FP_TYPE add = blob->energy + log(scale);
        if (nnos > 0) {
            for (int i = 0; i < nnos; ++i)
                frame->noise->spec[i] += add;
        }
        if (mb > 0.0f && freqwarp[0] < 3000.0f) {
            int idx = 0;
            do
                ++idx;
            while (freqwarp[idx] < 3000.0f);
            FP_TYPE rest = 1.0f - mb * 0.01f;
            for (int i = 0; i < idx; ++i)
                frame->noise->spec[i] += log_2(
                    (FP_TYPE)i * (1.0f / (FP_TYPE)idx * (1.0f - rest))
                    + rest + 1.0e-10f);
        }
    }

    free(faxis);
    free(amp);
    free(vt_freq);
    free(har_freq);
    free(lf_phse);
    free(lf_src);
    free(lf_dst);
    free(warped);

    if (args->Md != 0.0f) {
        FP_TYPE md = args->Md;
        FP_TYPE scale = exp_2(0.025f * md);
        FP_TYPE md01 = md * 0.01f;
        int nch = frame->noise->nchannel;

        for (int ch = 0; ch < nch; ++ch) {
            llsm_sinframe *eenv = frame->noise->eenv[ch];
            int nh = eenv->nhar;
            if (nh > 0) {
                for (int j = 0; j < nh; ++j)
                    eenv->ampl[j] *= scale;
            }
            frame->noise->emin[ch] *= 1.0f / scale;
        }
        if (md01 < 0.0f && nhar_new > 12) {
            for (int i = 12; i != nhar_new; ++i) {
                FP_TYPE fade = i > 30 ? 1.0f : (FP_TYPE)(i - 12) * 0.055555556f;
                FP_TYPE r = (FP_TYPE)rand() * 0.000030518509f - 0.5f;
                frame->vs_har_phse[i] += r * (fade * md01 * (FP_TYPE)(2.0 * M_PI));
            }
        }
    }
    return blob;
}

static void *lmix_shift_spectrum(
    lmix_frame *frame, FP_TYPE *freqwarp, resampler_args *args, FP_TYPE f0)
{
    return lmix_shift_spectrum_core(frame, freqwarp, args, f0, NULL);
}

static void *lmix_shift_spectrum_edge(
    lmix_frame *frame, FP_TYPE *freqwarp, resampler_args *args,
    FP_TYPE f0, void *anchor)
{
    return lmix_shift_spectrum_core(
        frame, freqwarp, args, f0, (lmix_shift_blob *)anchor);
}

static void lmix_grow_harmonics(lmix_frame *frame, int nhar)
{
    if (nhar <= frame->nhar)
        return;
    frame->vs_har_phse = (FP_TYPE *)realloc(
        frame->vs_har_phse, (size_t)nhar * sizeof(FP_TYPE));
    frame->vt_magn = (FP_TYPE *)realloc(
        frame->vt_magn,
        (size_t)nhar * (size_t)frame->vt_stride * sizeof(FP_TYPE));
}

/* The PE interpolates these phases in the frame loops.  Keep this expression
 * local rather than introducing the external lerp_phase float-call boundary. */
static inline FP_TYPE lmix_interp_phase(FP_TYPE a, FP_TYPE b, FP_TYPE t)
{
    return (FP_TYPE)atan2(
        (double)sin_2(b) + ((double)sin_2(a) - sin_2(b)) * t,
        (double)cos_2(b) + ((double)cos_2(a) - cos_2(b)) * t);
}

static void lmix_interp_frame(
    lmix_frame *dst, lmix_frame *a, lmix_frame *b, FP_TYPE t)
{
    int nhar_a = a->nhar;
    int nhar_b = b->nhar;
    int nhar = nhar_a >= nhar_b ? nhar_a : nhar_b;
    int nhar_min = nhar_a <= nhar_b ? nhar_a : nhar_b;

    if (a->f0 > 0.0f && b->f0 > 0.0f) {
        lmix_grow_harmonics(dst, nhar);
        dst->growl = b->growl + (a->growl - b->growl) * t;
        dst->distortion = b->distortion + (a->distortion - b->distortion) * t;
        dst->f0 = b->f0 + (a->f0 - b->f0) * t;
        dst->rd = b->rd + (a->rd - b->rd) * t;
        dst->vt_freq_step = b->vt_freq_step + (a->vt_freq_step - b->vt_freq_step) * t;
        if (nhar_min > 0) {
            for (int i = 0; i < nhar_min; ++i)
                dst->vs_har_phse[i] = lmix_interp_phase(
                    a->vs_har_phse[i], b->vs_har_phse[i], t);
        }
        if (nhar > nhar_min && nhar_a <= nhar_b) {
            for (int i = nhar_min; i < nhar; ++i)
                dst->vs_har_phse[i] = b->vs_har_phse[i];
        }
        int n = dst->vt_stride * nhar_min;
        if (n > 0) {
            for (int i = 0; i < n; ++i)
                dst->vt_magn[i] = b->vt_magn[i] +
                    (a->vt_magn[i] - b->vt_magn[i]) * t;
        }
        n = dst->vt_stride * nhar;
        if (n > dst->vt_stride * nhar_min && nhar_a <= nhar_b) {
            for (int i = dst->vt_stride * nhar_min; i < n; ++i)
                dst->vt_magn[i] = b->vt_magn[i];
        }
    } else if (a->f0 > 0.0f || b->f0 > 0.0f) {
        lmix_frame *src = a->f0 > 0.0f ? a : b;
        lmix_grow_harmonics(dst, nhar);
        dst->f0 = src->f0;
        dst->rd = src->rd;
        dst->vt_freq_step = src->vt_freq_step;
        dst->growl = src->growl;
        dst->distortion = src->distortion;
        if (src == b) {
            dst->vt_stride = b->vt_stride;
            dst->vs_har_phse = (FP_TYPE *)realloc(
                dst->vs_har_phse, (size_t)nhar * sizeof(FP_TYPE));
            dst->vt_magn = (FP_TYPE *)realloc(
                dst->vt_magn,
                (size_t)dst->vt_stride * (size_t)nhar * sizeof(FP_TYPE));
            if (nhar_b > 0) {
                for (int i = 0; i < nhar_b; ++i)
                    dst->vs_har_phse[i] = b->vs_har_phse[i];
            }
            int n = dst->vt_stride * nhar_b;
            for (int i = 0; i < n; ++i)
                dst->vt_magn[i] = b->vt_magn[i];
        }
    } else {
        dst->f0 = 0.0f;
        dst->rd = 1.0f;
        dst->vt_freq_step = 0.0f;
        dst->growl = 0.0f;
        dst->distortion = 0.0f;
        nhar = nhar_a >= nhar_b ? nhar_a : nhar_b;
    }
    dst->nhar = nhar;

    llsm_nosframe *nd = dst->noise;
    llsm_nosframe *nb = b->noise;
    if (nd->nnos > 0) {
        for (int i = 0; i < nd->nnos; ++i)
            nd->spec[i] = nb->spec[i] + (nd->spec[i] - nb->spec[i]) * t;
    }
    if (nd->nchannel > 0) {
        for (int ch = 0; ch < nd->nchannel; ++ch) {
            /* dst already contains a: blend the right frame b into that copy. */
            llsm_sinframe *ea = nb->eenv[ch];
            llsm_sinframe *eb = nd->eenv[ch];
            nd->emin[ch] = nb->emin[ch] + (nd->emin[ch] - nb->emin[ch]) * t;
            int nh = ea->nhar;
            if (eb->nhar < nh)
                nh = eb->nhar;
            if (nh <= 0)
                continue;
            for (int i = 0; i < nh; ++i) {
                eb->ampl[i] = ea->ampl[i] + (eb->ampl[i] - ea->ampl[i]) * t;
                eb->phse[i] = lmix_interp_phase(eb->phse[i], ea->phse[i], t);
            }
        }
    }
    dst->wave = NULL;
}

/* Log-frequency endpoints are single-precision model parameters. */
lmix_warp_bounds lmix_log_frequency_bounds(
    FP_TYPE low, FP_TYPE high, FP_TYPE warp)
{
    lmix_warp_bounds bounds;
    bounds.low = 5000.0 * log(1.0 + low / warp);
    bounds.high = 5000.0 * log(1.0 + high / warp);
    return bounds;
}

lmix_data *lmix_modify_data(
    lmix_data *lmix, resampler_args *args, const lmix_time_mapping mapping)
{
    // The input timing controls stay immutable; nmap tracks the working map.
    int nmap = mapping.nmap;
    int fs = g_output_sampling_rate;
    int hop = hop_from_fs(fs);
    FP_TYPE thop_out = (FP_TYPE)hop / (FP_TYPE)fs;
    int nfrm_out = (int)((double)fs * (double)mapping.length / (double)hop);
    int nfrm_in = lmix->conf.nfrm;
    int nmap1 = 0;
    FP_TYPE *src_map = NULL;
    FP_TYPE *dst_map = NULL;

    if (nmap == 0) {
        FP_TYPE left = mapping.offset + mapping.fixed_length;
        FP_TYPE right = mapping.total_dur - mapping.stretch_end;
        FP_TYPE leftover = right - left;
        FP_TYPE stretched = mapping.length - mapping.stretch;
        int ext_method;
        if (args->e == 1)
            ext_method = 1;
        else if (args->e == -1)
            ext_method = 2;
        else if (g_synthesis_duration_extension_method != 0)
            ext_method = g_synthesis_duration_extension_method;
        else if (stretched > 6.0f * leftover && leftover < 0.3f)
            ext_method = 1;
        else if (stretched < 0.0f)
            ext_method = 1;
        else
            ext_method = 2;

        log_file_printf("generate_time_mapping: ext_method = %d\n", ext_method);
        if (ext_method == 1) {
            src_map = (FP_TYPE *)calloc(3, sizeof(FP_TYPE));
            dst_map = (FP_TYPE *)calloc(3, sizeof(FP_TYPE));
            src_map[0] = mapping.offset;
            src_map[1] = left;
            src_map[2] = right;
            dst_map[0] = 0.0f;
            dst_map[1] = mapping.stretch;
            dst_map[2] = mapping.length;
            nmap = 3;
            nmap1 = 4;
            src_map = (FP_TYPE *)realloc(src_map, 4 * sizeof(FP_TYPE));
            dst_map = (FP_TYPE *)realloc(dst_map, 4 * sizeof(FP_TYPE));
        } else {
            int nrep = (int)ceil((double)stretched / leftover);
            int n = nrep + 2;
            src_map = (FP_TYPE *)calloc((size_t)n, sizeof(FP_TYPE));
            dst_map = (FP_TYPE *)calloc((size_t)n, sizeof(FP_TYPE));
            src_map[0] = mapping.offset;
            for (int i = 1; i < n; ++i)
                src_map[i] = (i & 1) ? left : right;
            dst_map[0] = 0.0f;
            dst_map[1] = mapping.stretch;
            for (int i = 2; i < n; ++i)
                dst_map[i] = dst_map[i - 1] + leftover;
            nmap = n;
            nmap1 = nrep + 3;
            src_map = (FP_TYPE *)realloc(
                src_map, (size_t)nmap1 * sizeof(FP_TYPE));
            dst_map = (FP_TYPE *)realloc(
                dst_map, (size_t)nmap1 * sizeof(FP_TYPE));
        }
        memmove(src_map + 1, src_map, (size_t)nmap * sizeof(FP_TYPE));
        memmove(dst_map + 1, dst_map, (size_t)nmap * sizeof(FP_TYPE));
    } else {
        src_map = (FP_TYPE *)calloc((size_t)nmap, sizeof(FP_TYPE));
        dst_map = (FP_TYPE *)calloc((size_t)nmap, sizeof(FP_TYPE));
        if (nmap > 0) {
            memcpy(src_map, mapping.map_src, (size_t)nmap * sizeof(FP_TYPE));
            memcpy(dst_map, mapping.map_dst, (size_t)nmap * sizeof(FP_TYPE));
        }
        nmap1 = nmap + 1;
        src_map = (FP_TYPE *)realloc(
            src_map, (size_t)nmap1 * sizeof(FP_TYPE));
        dst_map = (FP_TYPE *)realloc(
            dst_map, (size_t)nmap1 * sizeof(FP_TYPE));
        if (nmap > 0) {
            memmove(src_map + 1, src_map, (size_t)nmap * sizeof(FP_TYPE));
            memmove(dst_map + 1, dst_map, (size_t)nmap * sizeof(FP_TYPE));
        }
    }
    src_map[0] = src_map[1] - 1.0f;
    dst_map[0] = dst_map[1] - 1.0f;

    FP_TYPE *src_time = (FP_TYPE *)calloc((size_t)nfrm_in, sizeof(FP_TYPE));
    FP_TYPE *src_idx = (FP_TYPE *)malloc((size_t)nfrm_in * sizeof(FP_TYPE));
    if (nfrm_in > 0) {
        for (int i = 0; i < nfrm_in; ++i)
            src_idx[i] = (FP_TYPE)i;
        for (int i = 0; i < nfrm_in; ++i)
            src_time[i] = (FP_TYPE)i * lmix->conf.thop;
    }
    FP_TYPE *out_time = (FP_TYPE *)calloc((size_t)nfrm_out, sizeof(FP_TYPE));
    for (int i = 0; i < nfrm_out; ++i)
        out_time[i] = (FP_TYPE)i * thop_out;

    log_file_printf("modify_lmix_data: time mapping:\n");
    for (int i = 0; i < nmap1; ++i)
        log_file_printf("%7.3f -> %7.3f\n", src_map[i], dst_map[i]);

    FP_TYPE *src_at_out = cig_interp(dst_map, src_map, nmap1, out_time, nfrm_out);
    FP_TYPE *frm_at_out = cig_interp(
        src_time, src_idx, nfrm_in, src_at_out, nfrm_out);

    llsm_conf conf = lmix->conf;
    conf.nfrm = nfrm_out;
    conf.nhop = hop;
    conf.thop = thop_out;
    lmix_data *out = lmix_create(conf);
    log_file_printf("modify_lmix_data: interpolating lmix frames...\n");

    {
        int prev = nfrm_out > 0 ? (int)frm_at_out[0] : 0;
        for (int i = 0; i < nfrm_out; ++i) {
            FP_TYPE pos = frm_at_out[i];
            if (pos < 0.0f)
                pos = 0.0f;
            if (pos >= (FP_TYPE)(nfrm_in - 1))
                pos = (FP_TYPE)nfrm_in - 1.01f;
            int ia = (int)pos;
            lmix_frame *fa = lmix->frames[ia];
            lmix_frame *fb = lmix->frames[ia + 1];
            out->frames[i] = lmix_create_frame(
                out->conf, g_lmix_nfft, fb->nhar);
            lmix_copy_frame(out->frames[i], fa);
            if (ia == prev) {
                FP_TYPE frac = 1.0f - (FP_TYPE)fmod((double)pos, 1.0);
                lmix_interp_frame(out->frames[i], fa, fb, frac);
            }
            prev = ia;
        }
    }

    free(src_map);
    free(dst_map);
    free(src_at_out);
    free(frm_at_out);
    free(src_time);
    free(src_idx);

    log_file_printf("modify_lmix_data: generating f0 contour...\n");
    FP_TYPE pb_step = args->unknow12;
    if (pb_step <= 0.0f)
        pb_step = 0.625f / args->tempo;
    int npb = args->pitch_bend_size;
    FP_TYPE *pb_idx = (FP_TYPE *)malloc((size_t)npb * sizeof(FP_TYPE));
    FP_TYPE *pb_time = (FP_TYPE *)calloc((size_t)npb, sizeof(FP_TYPE));
    for (int i = 0; i < npb; ++i) {
        pb_idx[i] = (FP_TYPE)i;
        pb_time[i] = (FP_TYPE)i * pb_step;
    }
    FP_TYPE *pb_at_out = cig_interp(pb_time, pb_idx, npb, out_time, nfrm_out);
    FP_TYPE *pb_val = cig_interp(
        pb_idx, args->pitch_bend, npb, pb_at_out, nfrm_out);

    FP_TYPE *f0_src = (FP_TYPE *)calloc((size_t)nfrm_in, sizeof(FP_TYPE));
    FP_TYPE *f0_out = (FP_TYPE *)calloc((size_t)nfrm_out, sizeof(FP_TYPE));
    FP_TYPE *f0_note = (FP_TYPE *)calloc((size_t)nfrm_out, sizeof(FP_TYPE));
    for (int i = 0; i < nfrm_in; ++i) {
        if (lmix->frames[i] != NULL)
            f0_src[i] = lmix->frames[i]->f0;
    }
    for (int i = 0; i < nfrm_out; ++i)
        f0_out[i] = out->frames[i]->f0;

    FP_TYPE note_hz = midi_to_hz(midi_from_note(args->pitch));
    FP_TYPE t_ratio = (FP_TYPE)pow(2.0, 0.0008333333333333334 * args->t);
    note_hz *= t_ratio;

    FP_TYPE *voiced = (FP_TYPE *)calloc((size_t)nfrm_out + 2, sizeof(FP_TYPE));
    FP_TYPE *voiced_t = (FP_TYPE *)calloc((size_t)nfrm_out + 2, sizeof(FP_TYPE));
    int nvoiced = 1;
    for (int i = 0; i < nfrm_out; ++i) {
        if (f0_out[i] != 0.0f) {
            voiced[nvoiced] = f0_out[i];
            voiced_t[nvoiced] = (FP_TYPE)i;
            ++nvoiced;
        }
    }
    voiced[0] = voiced[1];
    voiced[nvoiced] = voiced[nvoiced - 1];
    voiced_t[0] = 0.0f;
    voiced_t[nvoiced] = (FP_TYPE)nfrm_out;
    nvoiced += 2;

    FP_TYPE *unit = (FP_TYPE *)malloc((size_t)nfrm_out * sizeof(FP_TYPE));
    for (int i = 0; i < nfrm_out; ++i)
        unit[i] = (FP_TYPE)i;
    FP_TYPE *tmp = cig_interp(voiced_t, voiced, nvoiced, unit, nfrm_out);
    free(unit);
    free(voiced);
    free(voiced_t);
    FP_TYPE *f0_smooth = cig_moving_avg(tmp, nfrm_out, 0.1f / out->conf.thop);

    FP_TYPE *noise = (FP_TYPE *)calloc((size_t)nfrm_out, sizeof(FP_TYPE));
    for (int i = 0; i < nfrm_out; ++i)
        noise[i] = randn(0.0f, 1.0f) * args->MP * 0.01;
    FP_TYPE *noise_s = cig_moving_avg(noise, nfrm_out, 0.2f / out->conf.thop);

    for (int i = 0; i < nfrm_out; ++i) {
        FP_TYPE srcf0 = out->frames[i]->f0;
        FP_TYPE target;
        if (args->unknow12 > 0.0f && args->_s > 0.0f) {
            target = srcf0 * args->_s;
            f0_note[i] = target;
            f0_out[i] = target == 0.0f ? f0_smooth[i] : target;
            continue;
        }
        if (args->_s > 0.0f) {
            target = srcf0 * args->_s;
            f0_note[i] = target;
            f0_out[i] = target == 0.0f ? f0_smooth[i] : target;
            continue;
        }
        if (args->unknow12 > 0.0f)
            target = pb_val[i];
        else
            target = (FP_TYPE)pow(
                2.0, 0.0008333333333333334 * pb_val[i]) * note_hz;
        f0_note[i] = target;
        if (srcf0 <= 0.0f) {
            f0_out[i] = 0.0f;
        } else {
            FP_TYPE mix = target + (srcf0 - f0_smooth[i]) * args->modulation *
                (note_hz * 0.01f) / f0_smooth[i];
            f0_out[i] = exp(log(mix) + noise_s[i]);
        }
    }
    free(tmp);
    free(f0_smooth);
    free(noise);
    free(noise_s);
    free(pb_at_out);
    free(pb_val);
    free(out_time);
    free(pb_idx);
    free(pb_time);

    log_file_printf("modify_lmix_data: shifting pitch  ...\n");
    FP_TYPE *freqwarp = lmix_pitch_freqwarp(
        0.0f, out->conf.nosf, out->conf.nnos, out->conf.noswarp);
    int *runs = (int *)calloc((size_t)nfrm_out, 2 * sizeof(int));
    int nrun = 0;
    {
        int i = 0;
        while (i < nfrm_out) {
            if (f0_out[i] > 0.0f) {
                int end = i;
                while (end < nfrm_out && f0_out[end] > 0.0f)
                    ++end;
                runs[2 * nrun] = i;
                runs[2 * nrun + 1] = end;
                ++nrun;
                i = end + 1;
            } else {
                ++i;
            }
        }
    }
    for (int j = 0; j < nrun; ++j) {
        int begin = runs[2 * j];
        int end = runs[2 * j + 1];
        int span = end - begin;
        if (span > 3) {
            int mid0;
            int mid1;
            if (span <= 10) {
                mid0 = begin + span / 3;
                mid1 = end - span / 3;
            } else {
                mid0 = begin + 5;
                mid1 = end - 5;
            }
            void *left = lmix_shift_spectrum(
                out->frames[mid0], freqwarp, args, f0_out[mid0]);
            void *right = lmix_shift_spectrum(
                out->frames[mid1], freqwarp, args, f0_out[mid1]);
#pragma omp parallel for
            for (int k = mid0 + 1; k < mid1; ++k) {
                void *p = lmix_shift_spectrum(
                    out->frames[k], freqwarp, args, f0_out[k]);
                lmix_shift_blob_free(p);
            }
            for (int i = begin; i < mid0; ++i) {
                void *p = lmix_shift_spectrum_edge(
                    out->frames[i], freqwarp, args, f0_out[i], left);
                lmix_shift_blob_free(p);
            }
            for (int i = mid1 + 1; i < end; ++i) {
                void *p = lmix_shift_spectrum_edge(
                    out->frames[i], freqwarp, args, f0_out[i], right);
                lmix_shift_blob_free(p);
            }
            lmix_shift_blob_free(left);
            lmix_shift_blob_free(right);
        } else {
            for (int i = begin; i < end; ++i) {
                void *p = lmix_shift_spectrum(
                    out->frames[i], freqwarp, args, f0_out[i]);
                lmix_shift_blob_free(p);
            }
        }
    }
    {
        int max_nhar = out->conf.nhar;
        for (int i = 0; i < nfrm_out; ++i) {
            if (out->frames[i]->nhar > max_nhar)
                max_nhar = out->frames[i]->nhar;
        }
        out->conf.nhar = max_nhar;
    }

    if (args->A != 0.0f) {
        log_file_printf(
            "modify_lmix_data: applying 'A' flag (A = %5.2f)...\n",
            args->A);
        for (int i = 0; i < nfrm_out - 1; ++i) {
            FP_TYPE r = (FP_TYPE)pow(
                2.0,
                0.00001 * args->A *
                    ((hz_to_mel(f0_note[i + 1]) -
                      hz_to_mel(f0_note[i])) *
                     (1.0 / thop_out)));
            if (r > 1.2f)
                r = 1.2f;
            if (r < 0.5f)
                r = 0.5f;
            FP_TYPE gain = log_2(r + 1.0e-10f);
            gain *= 0.6931471805599453094f; // ln(r) = log2(r) * ln(2)
            lmix_frame *fr = out->frames[i];
            if (fr->f0 > 0.0f) {
                int n = fr->nhar * fr->vt_stride;
                for (int j = 0; j < n; ++j)
                    fr->vt_magn[j] += gain;
            }
            int n = fr->noise->nnos;
            for (int j = 0; j < n; ++j)
                fr->noise->spec[j] += gain;
        }
    }

    if (args->b != 0.0f) {
        log_file_printf(
            "modify_lmix_data: applying 'b' flag (b = %5.2f)...\n",
            args->b);
        FP_TYPE scale = 0.05f * args->b + 1.0f;
        FP_TYPE log_gain = log_2(scale + 1.0e-10f) * 0.6931471805599453094f;
        for (int i = 0; i < nfrm_out; ++i) {
            int span = sign_change_span(f0_out, nfrm_out, i);
            lmix_frame *fr = out->frames[i];
            if (!(f0_out[i] > 0.0f && span > 9) && fr->noise->nnos > 0) {
                if (span <= 4) {
                    for (int j = 0; j < fr->noise->nnos; ++j)
                        fr->noise->spec[j] += log_gain;
                } else {
                    FP_TYPE g;
                    if (f0_out[i] > 0.0f)
                        g = log_2((10.0f - (FP_TYPE)span) *
                            (args->b * 0.008333333333333333f) +
                            1.0000000001f) * 0.6931471805599453094f;
                    else
                        g = log_gain;
                    for (int j = 0; j < fr->noise->nnos; ++j)
                        fr->noise->spec[j] += g;
                }
            }
            if (fr->wave != NULL && fr->wave->nsample > 0) {
                for (int j = 0; j < fr->wave->nsample; ++j)
                    fr->wave->samples[j] *= scale;
            }
        }
    }

    if (args->Ms != 0) {
        FP_TYPE *phase = (FP_TYPE *)calloc((size_t)nfrm_out, sizeof(FP_TYPE));
        int nh = out->conf.nhar;
        if (nh > 20)
            nh = 20;
        for (int h = 0; h < nh; ++h) {
            for (int i = 0; i < nfrm_out; ++i) {
                lmix_frame *fr = out->frames[i];
                if (fr->nhar <= h) {
                    phase[i] = 0.0f;
                    continue;
                }
                double turns = (double)fr->vs_har_phse[h] * 0.1591549430918953;
                phase[i] = fr->vs_har_phse[h] -
                    (FP_TYPE)round(turns) * (FP_TYPE)(2.0 * M_PI);
            }
            lmix_unwrap_phase_track(
                phase, f0_out, nfrm_out, 30, 12 - args->Ms);
            for (int i = 0; i < nfrm_out; ++i) {
                if (out->frames[i]->nhar > h)
                    out->frames[i]->vs_har_phse[h] = phase[i];
            }
        }
        free(phase);
    }

    if (args->MC != 0.0f) {
        int nfft = g_lmix_nfft;
        int nbin = nfft / 2 + 1;
        FP_TYPE *faxis = llsm_uniform_faxis(nfft, fs);
        FP_TYPE x_hz[7] = {
            0.0f, 400.0f, 1000.0f, 2000.0f, 4000.0f, 8000.0f, 192000.0f
        };
        FP_TYPE y_db[7] = {
            0.0f, -5.0f, -2.5f, -0.5f, 0.0f, 0.0f, 0.0f
        };
        FP_TYPE *tilt = cig_interp(x_hz, y_db, 7, faxis, nbin);
        for (int i = 0; i < nbin; ++i)
            tilt[i] = exp_2(tilt[i]);

        for (int i = 0; i < nfrm_out; ++i) {
            lmix_frame *fr = out->frames[i];
            int nh = fr->nhar;
            if (nh <= 0)
                continue;
            FP_TYPE *amp = lmix_frame_get_harmonic_amp(fr, NULL);
            FP_TYPE *spec = cig_interpu(
                f0_out[i], (FP_TYPE)(nh + 1) * f0_out[i],
                amp, nh, faxis, nbin);
            for (int k = 0; k < nbin; ++k)
                spec[k] *= tilt[k];
            FP_TYPE *env = llsm_geometric_envelope(
                spec, nfft, fs, freqwarp, fr->noise->nnos);
            if (fr->noise->nnos > 0 &&
                freqwarp[0] <= f0_out[i] * (FP_TYPE)nh) {
                for (int k = 0; k < fr->noise->nnos; ++k) {
                    if (freqwarp[k] > f0_out[i] * (FP_TYPE)nh)
                        break;
                    FP_TYPE v = fr->noise->spec[k];
                    if (freqwarp[k] < 400.0f)
                        v -= 0.05f * args->MC;
                    FP_TYPE floorv = log_2(env[k] * args->MC * 0.005f) *
                        0.6931471805599453094f;
                    if (v < floorv)
                        v = floorv;
                    fr->noise->spec[k] = v;
                }
            }
            free(env);
            free(spec);
            free(amp);
        }
        free(tilt);
        free(faxis);
    }

    if (args->MG != 0.0f) {
        FP_TYPE v = args->MG * 0.01f;
        for (int i = 0; i < nfrm_out; ++i)
            out->frames[i]->growl = v;
    }
    if (args->MD != 0.0f) {
        FP_TYPE v = args->MD * 0.01f;
        for (int i = 0; i < nfrm_out; ++i)
            out->frames[i]->distortion = v;
    }

    free(runs);
    free(freqwarp);
    free(f0_src);
    free(f0_out);
    free(f0_note);
    return out;
}

void lmix_apply_envelope(lmix_data *lmix, const lmix_envelope *env)
{
    int nfrm = lmix->conf.nfrm;
    FP_TYPE thop = (FP_TYPE)lmix->conf.nhop / (FP_TYPE)g_output_sampling_rate;
    FP_TYPE *times = (FP_TYPE *)calloc((size_t)nfrm, sizeof(FP_TYPE));
    if (nfrm > 0) {
        for (int i = 0; i != nfrm; ++i)
            times[i] = (FP_TYPE)i * thop - env->t0;
    }
    FP_TYPE tx[7];
    FP_TYPE ty[7];

    tx[0] = 0.0f;
    tx[1] = env->p1;
    tx[2] = env->p1 + env->p2;
    tx[3] = env->p1 + env->p2 + env->p5;
    tx[6] = env->duration;
    tx[5] = env->duration - env->p4;
    tx[4] = tx[5] - env->p3;
    ty[0] = 0.0f;
    ty[1] = env->v1;
    ty[2] = env->v2;
    ty[3] = env->v5;
    ty[4] = env->v3;
    ty[5] = env->v4;
    ty[6] = 0.0f;

    FP_TYPE *gain = cig_interp(tx, ty, 7, times, nfrm);
    if (nfrm <= 0) {
        free(gain);
        free(times);
        return;
    }

    for (int i = 0; i != nfrm; ++i) {
        lmix_frame *fr = lmix->frames[i];
        FP_TYPE g = log_2(0.01f * gain[i] + 1.0e-10f);
        if (fr->f0 > 0.0f) {
            int n = fr->vt_stride * fr->nhar;
            if (n > 0) {
                for (int k = 0; k < n; ++k)
                    fr->vt_magn[k] += g;
            }
        }
        int n = fr->noise->nnos;
        if (n > 0) {
            for (int k = 0; k < n; ++k)
                fr->noise->spec[k] += g;
        }
    }
    free(gain);
    free(times);

    for (int i = 0; i != nfrm; ++i) {
        llsm_wave *wave = lmix->frames[i]->wave;
        if (wave == NULL)
            continue;
        int nsample = wave->nsample;
        FP_TYPE t0 = (FP_TYPE)i * thop - env->t0 + wave->t0;
        FP_TYPE *st = (FP_TYPE *)malloc((size_t)nsample * sizeof(FP_TYPE));
        if (nsample > 0) {
            FP_TYPE inv_fs = 1.0f / (FP_TYPE)wave->fs;
            for (int j = 0; j != nsample; ++j)
                st[j] = (FP_TYPE)j;
            for (int j = 0; j < nsample; ++j)
                st[j] = st[j] * inv_fs + t0;
        }
        FP_TYPE *sg = cig_interp(tx, ty, 7, st, nsample);
        if (nsample > 0) {
            for (int j = 0; j != nsample; ++j)
                wave->samples[j] = sg[j] * wave->samples[j] * 0.01f;
        }
        free(st);
        free(sg);
    }
}

void lmix_scale_volume(lmix_data *lmix, float scale)
{
    int nfrm = lmix->conf.nfrm;
    if (nfrm <= 0)
        return;
    FP_TYPE g = log_2(scale + 1.0e-10f);
    for (int i = 0; i < nfrm; ++i) {
        lmix_frame *fr = lmix->frames[i];
        if (g_synthesis_utau_style_normalization == 1) {
            FP_TYPE ng;
            if (fr->f0 <= 0.0f) {
                ng = 1.0000001e-10f;
            } else {
                ng = g;
                int n = fr->nhar * fr->vt_stride;
                if (n > 0) {
                    for (int k = 0; k < n; ++k)
                        fr->vt_magn[k] += g;
                }
            }
            int n = fr->noise->nnos;
            if (n > 0) {
                for (int k = 0; k < n; ++k)
                    fr->noise->spec[k] += ng;
            }
        } else {
            if (fr->f0 > 0.0f) {
                int n = fr->vt_stride * fr->nhar;
                if (n > 0) {
                    for (int k = 0; k < n; ++k)
                        fr->vt_magn[k] += g;
                }
            }
            int n = fr->noise->nnos;
            if (n > 0) {
                for (int k = 0; k < n; ++k)
                    fr->noise->spec[k] += g;
            }
            if (fr->wave != NULL && fr->wave->nsample > 0) {
                FP_TYPE *s = fr->wave->samples;
                int ns = fr->wave->nsample;
                for (int k = 0; k < ns; ++k)
                    s[k] *= scale;
            }
        }
    }
}

// cig_gensins that stops when the next frequency is < 50 Hz
static FP_TYPE *lmix_gensins_cutoff(
    FP_TYPE *freq, FP_TYPE *ampl, FP_TYPE *phse,
    int nsin, int fs, int n)
{
    int nuse = nsin;
    for (int i = 1; i < nsin; ++i) {
        if (freq[i] < 50.0f) {
            nuse = i;
            break;
        }
    }
    return cig_gensins(freq, ampl, phse, nuse, fs, n);
}

void lmix_normalize_data(lmix_data *lmix, float scale, float P)
{
    int nfrm = lmix->conf.nfrm;
    log_file_printf("normalize_lmix_data: pre-synthesizing...\n");
    int nslot = (nfrm + 1) / 2;
    int best = -1;
    int nvoiced = 0;
    FP_TYPE best_e = 0.0f;

#pragma omp parallel for
    for (int i = 0; i < nslot; ++i) {
        int fi = 2 * i;
        if (fi >= nfrm)
            continue;
        lmix_frame *fr = lmix->frames[fi];
        if (fr->nhar <= 0)
            continue;
        int nhar_save = fr->nhar;
        if (nhar_save != 1)
            fr->nhar = nhar_save >> 1;
        FP_TYPE *amp = lmix_frame_get_harmonic_amp(fr, NULL);
        FP_TYPE e = 0.0f;
        if (fr->nhar > 0) {
            for (int k = 0; k < fr->nhar; ++k)
                e += amp[k] * amp[k];
        }
        fr->nhar = nhar_save;
#pragma omp critical
        {
            nvoiced += 1;
            if (e > best_e) {
                best_e = e;
                best = fi;
            }
        }
        free(amp);
    }

    if (best < 0 || nvoiced <= 4)
        return;

    log_file_printf("normalize_lmix_data: synthesizing frame...\n");
    {
        lmix_frame *fr = lmix->frames[best];
        int nhar = fr->nhar;
        int fs = g_output_sampling_rate;
        FP_TYPE *freq = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
        FP_TYPE *phse = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
        FP_TYPE *amp = lmix_frame_get_harmonic_amp(fr, phse);
        if (nhar > 0) {
            for (int k = 0; k != nhar; ++k)
                freq[k] = ((FP_TYPE)k + 1.0f) * fr->f0;
        }
        FP_TYPE *x = lmix_gensins_cutoff(freq, amp, phse, nhar, fs, 2048);
        FP_TYPE peak = 0.0f;
        for (int k = 0; k < 2048; ++k) {
            FP_TYPE a = (FP_TYPE)fabs((double)x[k]);
            if (peak < a)
                peak = a;
        }
        free(x);
        free(freq);
        free(phse);
        free(amp);
        log_file_printf("normalize_lmix_data: normalizing...\n");
        if (peak >= 1.0e-10f)
            lmix_scale_volume(
                lmix, (scale / peak - 1.0f) * (0.01f * P) + 1.0f);
    }
}

// hop-synchronous Hann OLA; depth 0.1 = low band, 0.5 = residual with MG/MD
static FP_TYPE *lmix_ola_band(
    llsm_output *out, FP_TYPE *x, lmix_data *lmix, FP_TYPE depth)
{
    static const FP_TYPE cycle[6] = {
        -0.89830506f, -0.69491524f, -0.18644071f,
         0.62711865f,  0.83050841f,  0.32203388f
    };
    int ny = out->ny;
    int nhop = lmix->conf.nhop;
    int nfrm = lmix->conf.nfrm;
    int ninst = out->ninst;
    int *instants = out->instants;
    FP_TYPE *acc = (FP_TYPE *)calloc((size_t)ny, sizeof(FP_TYPE));
    FP_TYPE *weight = (FP_TYPE *)calloc((size_t)ny, sizeof(FP_TYPE));

    for (int i = 0; i < ninst; ++i) {
        int center = instants[i];
        int half;
        if (i == 0) {
            half = instants[1] - center; // always reads instants[1]
        } else if (i > 1 && i < ninst - 1) {
            int fwd = instants[i + 1] - center;
            int bwd = center - instants[i - 1];
            half = fwd <= bwd ? fwd : bwd;
        } else {
            half = center - instants[i - 1];
        }
        int winlen = 2 * half;
        FP_TYPE *win = (FP_TYPE *)malloc((size_t)winlen * sizeof(FP_TYPE));
        if (winlen > 0) {
            for (int t = 0; t != winlen; ++t)
                win[t] = (1.0f - (FP_TYPE)cos(
                    (double)t * (2.0 * M_PI / (double)(winlen - 1)))) * 0.5f;
        }

        int fi = nfrm - 1;
        if (nhop != 0 && fi > center / nhop)
            fi = center / nhop;
        lmix_frame *fr = lmix->frames[fi];

        int jitter;
        if (depth > 0.2f) {
            jitter = (int)(((double)rand() * 0.000030518509
                + cos((double)i * 1.047197551196598) - 0.5)
                * (double)half * (0.03 * (double)fr->growl));
        } else {
            jitter = 0;
        }

        FP_TYPE pulsed = depth * cycle[i % 6];
        FP_TYPE g0 = (((FP_TYPE)rand() * 0.000030518509f - 0.5f) * 0.3f + pulsed)
            * (1.5f * fr->growl) + 1.0f;
        FP_TYPE sign = (i & 1) != 0 ? -1.5f : 0.5f;
        FP_TYPE gain = ((FP_TYPE)rand() * 0.000030518509f + sign)
            * (depth * 1.5f * fr->distortion) + g0;

        int t0 = center - half;
        int t1 = winlen - half + center;
        for (int t = t0; t != t1; ++t) {
            if (t >= 0 && t < ny - 1 && t + jitter < ny - 1) {
                FP_TYPE w = win[t - t0];
                weight[t] += w;
                acc[t + jitter] += w * x[t] * gain;
            }
        }
        free(win);
    }

    if (ny > 0) {
        for (int i = 0; i < ny; ++i) {
            FP_TYPE fill = (weight[i] <= 1.0f) ? (1.0f - weight[i]) * x[i] : 0.0f;
            acc[i] += fill;
        }
    }
    free(weight);
    return acc;
}

// ln(e^a + e^b) with -100 floor and 1e-10 clamp
static FP_TYPE lmix_log_sum_exp(FP_TYPE a, FP_TYPE b)
{
    if (a < -100.0f)
        a = -100.0f;
    if (b < -100.0f)
        b = -100.0f;
    FP_TYPE s = exp_2(a) + exp_2(b);
    if (s < 1.0e-10f)
        return -23.02585f;
    return log_2(s);
}

// wrap (x + 3π) mod 2π into [-π, π)
static FP_TYPE lmix_wrap_phase_3pi(FP_TYPE x)
{
    double y = (double)x + 9.424777960769379;
    y = fmod(y, 6.283185307179586);
    return (FP_TYPE)(y - 3.141592653589793);
}

void lmix_superimpose_add(
    lmix_data *dst, lmix_data *src,
    int begin, int overlap, int ltruncate, int flag)
{
    log_file_printf(
        "lmix_superimpose_add: begin=%d, overlap=%d, ltruncate=%d\n",
        begin, overlap, ltruncate);

    FP_TYPE *f0s = (FP_TYPE *)calloc((size_t)overlap, sizeof(FP_TYPE));
    FP_TYPE *phases = (FP_TYPE *)calloc((size_t)overlap, sizeof(FP_TYPE));

    int src_i = 0;
    if (ltruncate >= 0)
        src_i = ltruncate;
    int src_nfrm = src->conf.nfrm;
    if (src_i < src_nfrm) {
        int flag_n2 = flag & ~2;
        int ov_i = src_i - ltruncate;
        int dst_i = src_i + begin - ltruncate;
        while (src_i < src_nfrm) {
            FP_TYPE w_dst = 1.0f;
            if (overlap > ov_i)
                w_dst = (FP_TYPE)ov_i / (FP_TYPE)overlap;
            if (dst_i >= dst->conf.nfrm)
                break;
            if (dst_i < 0)
                goto next_frame;

            lmix_frame *df = dst->frames[dst_i];
            lmix_frame *sf;
            if (0.3 / (double)dst->conf.thop > (double)overlap) {
                sf = src->frames[src_i];
                FP_TYPE rd_d = df->rd;
                FP_TYPE rd_s = sf->rd;
                if ((rd_d < 0.15f || rd_s < 0.15f) &&
                    fabs((double)rd_d - (double)rd_s) > 0.3) {
                    if (rd_d < rd_s)
                        rd_d = rd_s;
                    df->rd = rd_d;
                    sf->rd = rd_d;
                }
            } else {
                sf = src->frames[src_i];
            }

            FP_TYPE w_src = 1.0f - w_dst;
            int nhar_dst = df->nhar;
            int nhar_src = sf->nhar;
            int mix_voiced = 0;
            int nhar_max;
            int nhar_min;
            int n;

            if (df->f0 <= 0.0f) {
                if (sf->f0 <= 0.0f) {
                    df->f0 = 0.0f;
                    df->rd = 1.0f;
                    df->vt_freq_step = 0.0f;
                    nhar_max = nhar_src >= nhar_dst ? nhar_src : nhar_dst;
                    df->growl = 0.0f;
                    df->distortion = 0.0f;
                    if (sf != NULL)
                        goto keep_dst_harmonics;
                    goto copy_src_harmonics;
                }
                nhar_max = nhar_src >= nhar_dst ? nhar_src : nhar_dst;
                df->f0 = sf->f0;
                df->rd = sf->rd;
                df->vt_freq_step = sf->vt_freq_step;
                    df->growl = sf->growl;
                    df->distortion = sf->distortion;
                goto copy_src_harmonics;
            } else if (sf->f0 <= 0.0f) {
                nhar_max = nhar_src >= nhar_dst ? nhar_src : nhar_dst;
                goto keep_dst_harmonics;
            }

            {
                FP_TYPE rd0 = df->rd;
                FP_TYPE step0 = df->vt_freq_step;

                nhar_max = nhar_dst >= nhar_src ? nhar_dst : nhar_src;
                nhar_min = nhar_dst <= nhar_src ? nhar_dst : nhar_src;
                df->growl =
                    sf->growl + (df->growl - sf->growl) * w_src;
                df->distortion =
                    sf->distortion + (df->distortion - sf->distortion) * w_src;
                df->f0 = sf->f0 + (df->f0 - sf->f0) * w_src;
                df->rd = sf->rd + (rd0 - sf->rd) * w_src;
                df->vt_freq_step =
                    sf->vt_freq_step + (step0 - sf->vt_freq_step) * w_src;
                if (nhar_dst < nhar_max) {
                    df->vs_har_phse = (FP_TYPE *)realloc(
                        df->vs_har_phse,
                        (size_t)nhar_max * sizeof(FP_TYPE));
                    df->vt_magn = (FP_TYPE *)realloc(
                        df->vt_magn,
                        (size_t)nhar_max * (size_t)df->vt_stride *
                            sizeof(FP_TYPE));
                }
                if (nhar_min > 0) {
                    for (int i = 0; i < nhar_min; ++i)
                        df->vs_har_phse[i] = lerp_phase(
                            df->vs_har_phse[i], sf->vs_har_phse[i], w_src);
                }
                if (nhar_max > nhar_min && nhar_dst <= nhar_src) {
                    for (int i = nhar_min; i < nhar_max; ++i)
                        df->vs_har_phse[i] = sf->vs_har_phse[i];
                }
                n = df->vt_stride * nhar_min;
                if (flag_n2 != 0) {
                    if (n > 0) {
                        for (int i = 0; i < n; ++i)
                            df->vt_magn[i] = lmix_log_sum_exp(
                                sf->vt_magn[i], df->vt_magn[i]);
                    }
                    if (df->vt_stride * nhar_max <= n) {
                        mix_voiced = 0;
                        df->nhar = nhar_max;
                        goto mix_noise_max;
                    }
                } else {
                    if (n > 0) {
                        for (int i = 0; i < n; ++i)
                            df->vt_magn[i] =
                                sf->vt_magn[i] +
                                (df->vt_magn[i] - sf->vt_magn[i]) * w_src;
                    }
                    if (n >= df->vt_stride * nhar_max) {
                        df->nhar = nhar_max;
                        mix_voiced = flag == 0;
                        goto mix_noise_linear;
                    }
                }
                if (nhar_dst <= nhar_src) {
                    for (int i = n; i < df->vt_stride * nhar_max; ++i)
                        df->vt_magn[i] = sf->vt_magn[i];
                }
                mix_voiced = flag == 0;
                df->nhar = nhar_max;
                if (flag_n2 != 0)
                    goto mix_noise_max;
                goto mix_noise_linear;
            }

        keep_dst_harmonics:
            if (flag != 2)
                goto set_nhar_flag;
            n = df->vt_stride * nhar_dst;
            if (n > 0) {
                FP_TYPE g = w_src;
                FP_TYPE lg;

                if (g < 1.0e-10f)
                    g = 1.0e-10f;
                lg = log_2(g);
                for (int i = 0; i < n; ++i)
                    df->vt_magn[i] += lg;
            }
            goto set_nhar_unvoiced_mix;

        copy_src_harmonics:
            df->vt_stride = sf->vt_stride;
            df->vs_har_phse = (FP_TYPE *)realloc(
                df->vs_har_phse, (size_t)nhar_max * sizeof(FP_TYPE));
            df->vt_magn = (FP_TYPE *)realloc(
                df->vt_magn,
                (size_t)nhar_max * (size_t)df->vt_stride * sizeof(FP_TYPE));
            if (nhar_src > 0) {
                for (int i = 0; i < nhar_src; ++i)
                    df->vs_har_phse[i] = sf->vs_har_phse[i];
            }
            n = df->vt_stride * nhar_src;
            if (n > 0) {
                for (int i = 0; i < n; ++i)
                    df->vt_magn[i] = sf->vt_magn[i];
                if (flag != 2)
                    goto set_nhar_flag;
                if (w_src <= 0.9999999999) {
                    FP_TYPE lg = log_2(1.0f - w_src);
                    for (int i = 0; i < n; ++i)
                        df->vt_magn[i] += lg;
                } else {
                    for (int i = 0; i < n; ++i)
                        df->vt_magn[i] -= 23.02585092994046f;
                }
                goto set_nhar_unvoiced_mix;
            }
            if (flag != 2)
                goto set_nhar_flag;

        set_nhar_unvoiced_mix:
            mix_voiced = 0;
            df->nhar = nhar_max;
            goto mix_noise_linear;

        set_nhar_flag:
            mix_voiced = flag == 0;
            df->nhar = nhar_max;
            if (flag_n2 != 0)
                goto mix_noise_max;

        mix_noise_linear:
            if (df->noise->nnos > 0) {
                for (int i = 0; i < df->noise->nnos; ++i)
                    df->noise->spec[i] =
                        sf->noise->spec[i] +
                        (df->noise->spec[i] - sf->noise->spec[i]) * w_src;
            }
            goto mix_noise_channels;

        mix_noise_max:
            if (df->noise->nnos > 0) {
                for (int i = 0; i < df->noise->nnos; ++i) {
                    FP_TYPE v = df->noise->spec[i];
                    if (v < sf->noise->spec[i])
                        v = sf->noise->spec[i];
                    df->noise->spec[i] = v;
                }
            }

        mix_noise_channels:
            for (int ch = 0; ch < df->noise->nchannel; ++ch) {
                llsm_sinframe *ea = sf->noise->eenv[ch];
                llsm_sinframe *eb = df->noise->eenv[ch];
                if (mix_voiced || (df->f0 > 0.0f && sf->f0 > 0.0f)) {
                    df->noise->emin[ch] =
                        sf->noise->emin[ch] +
                        (df->noise->emin[ch] - sf->noise->emin[ch]) * w_src;
                } else {
                    FP_TYPE v = df->noise->emin[ch];
                    if (v < sf->noise->emin[ch])
                        v = sf->noise->emin[ch];
                    df->noise->emin[ch] = v;
                }
                int nh = ea->nhar;
                if (eb->nhar <= nh)
                    nh = eb->nhar;
                if (nh > 0) {
                    for (int i = 0; i < nh; ++i) {
                        eb->ampl[i] =
                            ea->ampl[i] +
                            (eb->ampl[i] - ea->ampl[i]) * w_src;
                        eb->phse[i] = lerp_phase(
                            eb->phse[i], ea->phse[i], w_src);
                    }
                }
                if ((unsigned)(flag - 1) <= 1u &&
                    ((df->f0 <= 0.0f &&
                      df->noise->emin[ch] >
                          sf->noise->emin[ch] + sf->noise->emin[ch]) ||
                     (sf->f0 <= 0.0f &&
                      sf->noise->emin[ch] >
                          df->noise->emin[ch] + df->noise->emin[ch])) &&
                    nh > 0) {
                    memset(eb->ampl, 0, (size_t)nh * sizeof(FP_TYPE));
                }
            }
            if (overlap > ov_i)
                f0s[ov_i] = dst->frames[dst_i]->f0;
            src_nfrm = src->conf.nfrm;

        next_frame:
            ++src_i;
            ++ov_i;
            ++dst_i;
        }
    }

    src_nfrm = src->conf.nfrm;
    if (src_nfrm > 0) {
        int di = begin - ltruncate;
        for (int i = 0; i < src_nfrm; ++i) {
            llsm_wave *sw = src->frames[i]->wave;
            if (sw != NULL) {
                int idx = 0;
                if (di > 0) {
                    idx = dst->conf.nfrm - 1;
                    if (di < dst->conf.nfrm)
                        idx = di;
                }
                lmix_frame *df = dst->frames[idx];
                if (df->wave == NULL) {
                    df->wave = llsm_create_wave(
                        sw->samples, sw->nsample, sw->fs);
                    df->wave->t0 = sw->t0;
                }
            }
            ++di;
        }
    }

    if (begin >= 0 && overlap > 2 &&
        overlap + begin < dst->conf.nfrm && dst->conf.nhar > 0) {
        int harm = 0;
        int hmax;
        do {
            int ii;
            for (ii = 0; ii < overlap; ++ii) {
                for (;;) {
                    lmix_frame *fr = dst->frames[begin + ii];
                    f0s[ii] = fr->f0;
                    if (fr->nhar > harm)
                        break;
                    phases[ii] = 0.0f;
                    ++ii;
                    if (ii >= overlap)
                        goto unwrap_call;
                }
                phases[ii] = lmix_wrap_phase_3pi(
                    dst->frames[begin + ii]->vs_har_phse[harm]);
            }
        unwrap_call:
            lmix_unwrap_phase_track(phases, f0s, overlap, 30, 5);
            for (ii = 0; ii < overlap; ++ii) {
                lmix_frame *fr = dst->frames[begin + ii];
                if (fr->nhar > harm)
                    fr->vs_har_phse[harm] = phases[ii];
            }
            ++harm;
            hmax = 25;
            if (dst->conf.nhar <= 25)
                hmax = dst->conf.nhar;
        } while (hmax > harm);
    }

    free(f0s);
    free(phases);
}

FP_TYPE *lmix_render_wave(llsm_output *out, lmix_data *lmix)
{
    FP_TYPE cutoff = 1500.0f / out->fs;
    FP_TYPE *fir = cig_winfir(32, cutoff, 0.0f, "lowpass", "hanning");
    FP_TYPE *filt = cig_convolution(out->y, fir, out->ny, 32);
    int ny = out->ny;
    FP_TYPE *residual = (FP_TYPE *)calloc((size_t)ny, sizeof(FP_TYPE));
    if (ny > 0) {
        for (int i = 0; i < ny; ++i)
            residual[i] = out->y[i] - filt[i + 16];
    }
    FP_TYPE *low = lmix_ola_band(out, filt + 16, lmix, 0.1f);
    FP_TYPE *high = lmix_ola_band(out, residual, lmix, 0.5f);
    if (ny > 0) {
        for (int i = 0; i < ny; ++i)
            low[i] += high[i];
    }
    free(filt);
    free(residual);
    free(high);
    free(fir);
    return low;
}

