#include "nebula.h"

#include "../external/fastapprox/fasttrig.h"
#include "../external/ciglet/ciglet.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


// IF/SNR: heterodyne each band, two centred complex FIRs, residual power, interpolate
static void get_wcoef(const char *name, float a[4])
{
    if (strcmp(name, "nuttall83") == 0) {
        a[0] = 0.338946f; a[1] = 0.481973f; a[2] = 0.161054f; a[3] = 0.018027f;
    } else if (strcmp(name, "nuttall98") == 0) {
        a[0] = 0.3635819f; a[1] = 0.4891775f; a[2] = 0.1365995f; a[3] = 0.0106411f;
    } else if (strcmp(name, "hanning") == 0) {
        a[0] = 0.5f; a[1] = 0.5f; a[2] = 0.0f; a[3] = 0.0f;
    } else if (strcmp(name, "hamming") == 0) {
        a[0] = 0.54f; a[1] = 0.46f; a[2] = 0.0f; a[3] = 0.0f;
    } else if (strcmp(name, "blackman") == 0) {
        a[0] = 0.42f; a[1] = 0.5f; a[2] = 0.08f; a[3] = 0.0f;
    } else if (strcmp(name, "blackman-harris") == 0) {
        a[0] = 0.4243801f; a[1] = 0.4973406f; a[2] = 0.0782793f; a[3] = 0.0f;
    } else {
        a[0] = 1.0f; a[1] = 0.0f; a[2] = 0.0f; a[3] = 0.0f;
    }
}

// cosine-sum window, normalized to unit sum
static float *make_window(const float aw[4], int nw, float omegaw)
{
    float *w = (float *)calloc((size_t)nw, sizeof(*w));

    for (int k = 0; k < 4 && aw[k] > 0.0f; ++k)
        for (int j = 0; j < nw; ++j)
            w[j] += aw[k] * fastcosfull((j - nw / 2) * (k * omegaw));

    if (nw > 0) {
        float sum = 0.0f;
        for (int i = 0; i != nw; ++i)
            sum += w[i];
        sum = 1.0f / sum;
        for (int i = 0; i != nw; ++i)
            w[i] *= sum;
    }
    return w;
}

static float neb_cabs(float re, float im)
{
    return sqrtf(re * re + im * im);
}

// interleaved real/imag recurrence
static float *make_oscillator(int n, float omega)
{
    float *osc = (float *)malloc((size_t)(2 * n) * sizeof(*osc));
    float c = cos(omega);
    float c2 = c + c;

    osc[0] = 1.0f;
    osc[1] = 0.0f;
    osc[2] = c;
    osc[3] = sin(omega);
    for (int i = 2; i < n; ++i) {
        osc[2 * i] = osc[2 * (i - 1)] * c2 - osc[2 * (i - 2)];
        osc[2 * i + 1] = osc[2 * (i - 1) + 1] * c2
                       - osc[2 * (i - 2) + 1];
    }
    return osc;
}

// complex overlap-add convolution; output length nx + nh - 1
// Every input block is transformed up-front and kept in a small plan record,
// the (real) kernel is transformed once, then each block is multiplied in
// place, inverse transformed and overlap-added.
typedef struct fftconv_plan {
    float **re;
    float **im;
    int nfft;
    int nblock;
    int step;
} fftconv_plan;

static void
fftconv_complex(const float *xr, const float *xi, int nx,
                const float *hr, int nh, float *yr, float *yi)
{
    int nfft = (int)pow(2.0, ceil(log2((double)nh)) + 2.0);
    int step = nfft - nh;
    fftconv_plan *plan = (fftconv_plan *)malloc(sizeof(*plan));
    int nblock = nx / step + 1;
    float **bre = (float **)calloc((size_t)nblock, sizeof(float *));
    float **bim;
    float *work;
    float *kre;
    float *kim;
    float *hre;
    float *him;
    const float *sxr = xr;
    const float *sxi = xi;
    int outlen;
    int offset;
    int remain;
    int b;
    int i;

    if (nblock == 0) {
        plan->re = bre;
        bim = (float **)calloc((size_t)0, sizeof(float *));
    } else {
        for (unsigned int j = 0; j != (unsigned int)nblock; ++j)
            bre[j] = (float *)calloc((size_t)nfft, sizeof(float));
        plan->re = bre;
        bim = (float **)calloc((size_t)nblock, sizeof(float *));
        for (unsigned int j = 0; j != (unsigned int)nblock; ++j)
            bim[j] = (float *)calloc((size_t)nfft, sizeof(float));
    }
    plan->im = bim;
    work = (float *)malloc((size_t)(8 * nfft) * sizeof(float));
    if (nblock > 0) {
        offset = 0;
        remain = nx;
        for (b = 0; b != nblock; ++b) {
            int take;
            for (i = 0; i < step && i + offset < nx; ++i) {
                bre[b][i] = sxr[i];
                bim[b][i] = xi != NULL ? sxi[i] : 0.0f;
            }
            take = i;
            float *pbr = bre[b];
            float *pbi = bim[b];
            if (nfft > take) {
                memset(&pbr[take], 0, (size_t)(4 * (nfft - take)));
                memset(&pbi[take], 0, (size_t)(4 * (nfft - take)));
            }
            cig_fft(pbr, pbi, pbr, pbi, nfft, work, -1.0f);
            offset += step;
            remain -= step;
            sxr += step;
            sxi += step;
        }
    }
    free(work);
    plan->nfft = nfft;
    plan->nblock = nblock;
    plan->step = step;
    outlen = nh + nx - 1;
    work = (float *)malloc((size_t)(8 * nfft) * sizeof(float));
    if (outlen > 0) {
        memset(yr, 0, (size_t)(4 * (nh + nx) - 4));
        if (yi != NULL)
            memset(yi, 0, (size_t)(4 * (nh + nx) - 4));
    }
    kre = (float *)calloc((size_t)nfft, sizeof(float));
    kim = (float *)calloc((size_t)nfft, sizeof(float));
    hre = (float *)calloc((size_t)nfft, sizeof(float));
    him = (float *)calloc((size_t)nfft, sizeof(float));
    if (nh > 0)
        memcpy(hre, hr, (size_t)(4 * nh));
    cig_fft(hre, him, kre, kim, nfft, work, -1.0f);
    if (nblock > 0) {
        int idx = 0;
        int limit = nfft;
        float *srcr = hre;
        float *srci = him;
        for (b = 0; b != nblock; ++b) {
            float *pbr = bre[b];
            float *pbi = bim[b];
            if (nfft != 0) {
                for (i = 0; i != nfft; ++i) {
                    hre[i] = pbr[i] * kre[i] - pbi[i] * kim[i];
                    him[i] = kim[i] * pbr[i] + kre[i] * pbi[i];
                }
                cig_fft(hre, him, hre, him, nfft, work, 1.0f);
                if (yi != NULL) {
                    for (i = idx; i != limit; ++i) {
                        if (outlen > i) {
                            yr[i] += srcr[i];
                            yi[i] += srci[i];
                        }
                    }
                } else {
                    for (i = idx; i != limit; ++i) {
                        if (outlen > i)
                            yr[i] += srcr[i];
                    }
                }
            } else {
                cig_fft(hre, him, hre, him, 0, work, 1.0f);
            }
            limit += step;
            idx += step;
            srcr -= step;
            srci -= step;
        }
        nblock = plan->nblock;
    }
    free(kre);
    free(kim);
    free(hre);
    free(him);
    free(work);
    if (nblock != 0) {
        for (i = 0; i != nblock; ++i)
            free(plan->re[i]);
        b = 0;
        do {
            free(plan->im[b++]);
        } while (b != i);
    }
    free(plan->re);
    free(plan->im);
    free(plan);
}

// after first decimation: phase-normalise, second FIR, residual power, pad 256 zeros
static int x2as_decimated(float *yr, float *yi, int nvalid,
                          const float aw[4], const float aws[4],
                          float bandwidth_fc, int decim, float **out)
{
    const int nwork = nvalid + 256;
    *out = NULL;
    if (nvalid <= 0 || bandwidth_fc <= 0.0f || decim <= 0)
        return 0;
    int nh = (int)ceil(4.0 / (decim * bandwidth_fc));
    if (nh <= 0)
        return 0;
    float omega = decim * bandwidth_fc * (0.5 * M_PI);
    float *w = make_window(aw, nh, omega);
    float *ws = make_window(aws, nh, omega);
    float *y2r = (float *)calloc((size_t)(nwork + nh - 1), sizeof(*y2r));
    float *y2i = (float *)calloc((size_t)(nwork + nh - 1), sizeof(*y2i));
    float *power = (float *)calloc((size_t)nwork, sizeof(*power));
    float *smooth = (float *)calloc((size_t)(nwork + nh - 1), sizeof(*smooth));
    if (!w || !ws || !y2r || !y2i || !power || !smooth) {
        free(w); free(ws); free(y2r); free(y2i); free(power); free(smooth);
        return 0;
    }

    for (int i = 0; i < nwork; ++i) {
        float mag = neb_cabs(yr[i], yi[i]);
        float den = mag + 1.0e-10f;
        yr[i] /= den;
        yi[i] /= den;
    }
    fftconv_complex(yr, yi, nwork, w, nh, y2r, y2i);
    for (int i = 0; i < nwork; ++i) {
        int src = i + nh / 2;
        // second unit-normalize before the phase residual
        float mag = neb_cabs(y2r[src], y2i[src]);
        float den = mag + 1.0e-10f;
        float pr = y2r[src] / den;
        float pi = y2i[src] / den;
        float dr = yr[i] - pr;
        float di = yi[i] - pi;
        power[i] = dr * dr + di * di;
    }
    fftconv_complex(power, NULL, nwork, ws, nh, smooth, NULL);

    *out = (float *)malloc((size_t)nvalid * sizeof(**out));
    if (*out) {
        for (int i = 0; i < nvalid; ++i)
            (*out)[i] = smooth[i + nh / 2];
    }
    free(w); free(ws); free(y2r); free(y2i); free(power); free(smooth);
    return *out != NULL;
}

// carrier = heterodyne; bandwidth_fc = window/decimation. SNR2 only changes carrier to 2*fc
static int extract_band_snr(const float *x, int nx, float carrier,
                            float bandwidth_fc,
                            const float aw[4], const float aws[4],
                            const int *naxis, int nfrm, float *dst)
{
    if (!x || !naxis || !dst || nx <= 0 || carrier <= 0.0f ||
        bandwidth_fc <= 0.0f)
        return 0;
    int nw = (int)ceil(4.0 / bandwidth_fc);
    int decim = (int)floor(0.5 / bandwidth_fc);
    if (nw <= 0 || decim <= 0)
        return 0;
    int ndec = nx / decim;
    if (ndec <= 0)
        return 0;

    float *w = make_window(aw, nw, 0.5 * M_PI * bandwidth_fc);
    float *osc = make_oscillator(nx, (float)(2.0 * M_PI * carrier));
    float *mix_r = (float *)malloc((size_t)nx * sizeof(*mix_r));
    float *mix_i = (float *)malloc((size_t)nx * sizeof(*mix_i));
    float *conv_r = (float *)calloc((size_t)(nx + nw - 1), sizeof(*conv_r));
    float *conv_i = (float *)calloc((size_t)(nx + nw - 1), sizeof(*conv_i));
    float *dec_r = (float *)calloc((size_t)(ndec + 256), sizeof(*dec_r));
    float *dec_i = (float *)calloc((size_t)(ndec + 256), sizeof(*dec_i));
    if (!w || !osc || !mix_r || !mix_i || !conv_r || !conv_i ||
        !dec_r || !dec_i) {
        free(w); free(osc); free(mix_r); free(mix_i);
        free(conv_r); free(conv_i); free(dec_r); free(dec_i);
        return 0;
    }

    for (int i = 0; i < nx; ++i) {
        mix_r[i] = x[i] * osc[2 * i];
        mix_i[i] = x[i] * osc[2 * i + 1];
    }
    fftconv_complex(mix_r, mix_i, nx, w, nw, conv_r, conv_i);
    for (int i = 0; i < ndec; ++i) {
        int src = i * decim + nw / 2;
        dec_r[i] = conv_r[src];
        dec_i[i] = conv_i[src];
    }
    float *as = NULL;
    if (!x2as_decimated(dec_r, dec_i, ndec, aw, aws,
                         bandwidth_fc, decim, &as)) {
        free(w); free(osc); free(mix_r); free(mix_i);
        free(conv_r); free(conv_i); free(dec_r); free(dec_i);
        return 0;
    }
    float *axis = (float *)malloc((size_t)ndec * sizeof(*axis));
    float *query = (float *)malloc((size_t)nfrm * sizeof(*query));
    if (!axis || !query) {
        free(axis); free(query); free(as); free(w); free(osc);
        free(mix_r); free(mix_i); free(conv_r); free(conv_i); free(dec_r); free(dec_i);
        return 0;
    }
    for (int i = 0; i < ndec; ++i)
        axis[i] = (float)(i * decim);
    for (int i = 0; i < nfrm; ++i)
        query[i] = (float)naxis[i];
    float *interp = cig_interp(axis, as, ndec, query, nfrm);
    if (interp) {
        memcpy(dst, interp, (size_t)nfrm * sizeof(*dst));
        free(interp);
    }

    free(axis); free(query); free(as); free(w); free(osc);
    free(mix_r); free(mix_i); free(conv_r); free(conv_i); free(dec_r); free(dec_i);
    return interp != NULL;
}

// ciglet IF detector, clamped to [lower, upper]
static void estimate_if(const float *x, int nx, const int *naxis, int nfrm,
                        float fc, float fres, float lower, float upper, float *dst)
{
    if (!dst)
        return;
    ifdetector *ifd = cig_create_ifdetector(fc, fres);
    if (!ifd)
        return;
    for (int t = 0; t < nfrm; ++t) {
        int nh = ifd->nh;
        int begin = naxis[t] - nh / 2;
        float *frame = (float *)malloc((size_t)nh * sizeof(*frame));
        float value = fc;
        if (frame) {
            for (int i = 0; i < nh; ++i) {
                int src = begin + i;
                frame[i] = (src >= 0 && src < nx) ? x[src] : 0.0f;
            }
            value = cig_ifdetector_estimate(ifd, frame, nh);
            free(frame);
        }
        if (value < lower)
            value = lower;
        if (value > upper)
            value = upper;
        dst[t] = value;
    }
    cig_delete_ifdetector(ifd);
}

// per-band SNR/IF; SNR2 reuses the fc-sized window with carrier 2*fc
void estimate_ifsnr(float *x, int nx, const int *naxis, int nfrm,
                    const float *fc, const float *fres, int nband, const char *wtype,
                    const char *stype, float **SNR1, float **IF1,
                    float **SNR2, float **IF2)
{
    if (!x || !naxis || !fc || !fres || !SNR1 || !IF1 ||
        nx <= 0 || nfrm <= 0 || nband <= 0)
        return;
    float aw[4];
    float aws[4];
    get_wcoef(wtype ? wtype : "nuttall98", aw);
    get_wcoef(stype ? stype : "hanning", aws);

#pragma omp parallel for
    for (int b = 0; b < nband; ++b) {
        float f = fc[b];
        float bw = fres[b];
        float *snr1 = (float *)calloc((size_t)nfrm, sizeof(*snr1));
        float *snr2 = (SNR2 || IF2)
                          ? (float *)calloc((size_t)nfrm, sizeof(*snr2))
                          : NULL;
        if (snr1 && (!SNR2 && !IF2 ? 1 : snr2 != NULL)) {
            extract_band_snr(x, nx, f, f, aw, aws, naxis, nfrm, snr1);
            for (int t = 0; t < nfrm; ++t)
                SNR1[t][b] = snr1[t];

            if (SNR2) {
                extract_band_snr(x, nx, f * 2.0f, f, aw, aws,
                                 naxis, nfrm, snr2);
                for (int t = 0; t < nfrm; ++t)
                    SNR2[t][b] = snr2[t];
            }

            estimate_if(x, nx, naxis, nfrm, f, bw, f * 0.5f, f * 1.5f, snr1);
            for (int t = 0; t < nfrm; ++t)
                IF1[t][b] = snr1[t];

            if (IF2) {
                estimate_if(x, nx, naxis, nfrm, f * 2.0f, bw,
                            f * 1.5f, f * 2.5f, snr2);
                for (int t = 0; t < nfrm; ++t)
                    IF2[t][b] = snr2[t];
            }
        }
        free(snr1);
        free(snr2);
    }

}
