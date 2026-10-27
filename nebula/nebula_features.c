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

// after first decimation: phase-normalise, second FIR, residual power.
// In-place: the smoothed power envelope is written back into yr.
// Windows are made by the caller (analysis h, smoothing hs) and shared
// between both passes.
static void x2as_decimated(float *yr, float *yi, int nwork,
                           const float *h, const float *hs, int nh)
{
    int nh2 = nh / 2;
    float *mag = (float *)malloc((size_t)(4 * nwork));
    float *y2r, *y2i, *pw;
    int i;

    for (i = 0; i != nwork; ++i)
        mag[i] = neb_cabs(yr[i], yi[i]);
    for (i = 0; i != nwork; ++i) {
        yr[i] = yr[i] / (mag[i] + 1.0e-10);
        yi[i] = yi[i] / (mag[i] + 1.0e-10);
    }
    free(mag);

    y2r = (float *)calloc((size_t)(nh + nwork), sizeof(*y2r));
    y2i = (float *)calloc((size_t)(nh + nwork), sizeof(*y2i));
    fftconv_complex(yr, yi, nwork, h, nh, y2r, y2i);

    pw = (float *)malloc((size_t)(4 * nwork));
    for (i = 0; i != nwork; ++i)
        pw[i] = neb_cabs(y2r[nh2 + i], y2i[nh2 + i]);
    for (i = 0; i != nwork; ++i) {
        float pr = y2r[nh2 + i] / (pw[i] + 1.0e-10);
        float pi = y2i[nh2 + i] / (pw[i] + 1.0e-10);
        y2r[nh2 + i] = pr;
        y2i[nh2 + i] = pi;
        yr[i] = yr[i] - pr;
        yi[i] = yi[i] - pi;
        pw[i] = yr[i] * yr[i] + yi[i] * yi[i];
    }
    free(y2i);

    fftconv_complex(pw, NULL, nwork, hs, nh, y2r, NULL);

    memcpy(yr, y2r + nh2, (size_t)(4 * nwork));
    free(y2r);
    free(pw);
}


// per-band SNR/IF worker: heterodyne, decimate, x2as, interpolate, scatter.
// The SNR2 pass re-runs only the carrier at 2*fc; the three windows, the mix
// and decimation buffers are shared between both passes.
void estimate_ifsnr(float *x, int nx, const int *naxis, int nfrm,
                    const float *fc, const float *fres, int nband,
                    const char *wtype, const char *stype,
                    float **SNR1, float **IF1, float **SNR2, float **IF2)
{
    float aw[4];
    float aws[4];
    float *query = (float *)calloc((size_t)nfrm, sizeof(*query));
    for (int t = 0; t < nfrm; ++t)
        query[t] = (float)naxis[t];
    get_wcoef(wtype, aw);
    get_wcoef(stype, aws);

#pragma omp parallel for
    for (int b = 0; b < nband; ++b) {
        float f = fc[b];
        int decim = (int)floor(0.5 / f);
        int ndec = nx / decim;
        float *axis = (float *)calloc((size_t)ndec, sizeof(*axis));
        for (int i = 0; i < ndec; ++i)
            axis[i] = (float)(i * decim);

        float carrier = f * 6.283185307179586;
        float bw = fres[b];
        int nw = (int)ceil(4.0 / bw);
        double scaled_bw = decim * (double)bw;
        int nh = (int)ceil(4.0 / scaled_bw);
        float omega = scaled_bw * (0.5 * M_PI);
        float *w = make_window(aw, nw, 0.5 * M_PI * bw);
        float *wf = make_window(aw, nh, omega);
        float *ws = make_window(aws, nh, omega);
        float *mix_r = (float *)calloc((size_t)(nw + nx), sizeof(*mix_r));
        float *mix_i = (float *)calloc((size_t)(nw + nx), sizeof(*mix_i));
        float *conv_r, *conv_i;
        float *dec_r, *dec_i;
        float *osc = make_oscillator(nx, (float)carrier);
        int i, t;

        for (i = 0; i < nx; ++i) {
            mix_r[i] = osc[2 * i] * x[i];
            mix_i[i] = x[i] * osc[2 * i + 1];
        }
        conv_r = (float *)calloc((size_t)(nw + nx), sizeof(*conv_r));
        conv_i = (float *)calloc((size_t)(nw + nx), sizeof(*conv_i));
        fftconv_complex(mix_r, mix_i, nx, w, nw, conv_r, conv_i);
        if (nx > 0) {
            memcpy(mix_r, conv_r + nw / 2, (size_t)(4 * nx));
            memcpy(mix_i, conv_i + nw / 2, (size_t)(4 * nx));
        }
        free(conv_r);
        free(conv_i);
        free(osc);
        int nwork = ndec + 256;
        dec_r = (float *)calloc((size_t)nwork, sizeof(*dec_r));
        dec_i = (float *)calloc((size_t)nwork, sizeof(*dec_i));
        for (i = 0; i < ndec; ++i) {
            dec_r[i] = mix_r[i * decim];
            dec_i[i] = mix_i[i * decim];
        }
        x2as_decimated(dec_r, dec_i, nwork, wf, ws, nh);
        {
            float *interp = cig_interp(axis, dec_r, ndec, query, nfrm);
            for (t = 0; t < nfrm; ++t)
                SNR1[t][b] = interp[t];
            free(interp);
        }

        if (SNR2) {
            float *osc2 = make_oscillator(nx, (float)(carrier + carrier));
            for (i = 0; i < nx; ++i) {
                mix_r[i] = osc2[2 * i] * x[i];
                mix_i[i] = x[i] * osc2[2 * i + 1];
            }
            conv_r = (float *)calloc((size_t)(nw + nx), sizeof(*conv_r));
            conv_i = (float *)calloc((size_t)(nw + nx), sizeof(*conv_i));
            fftconv_complex(mix_r, mix_i, nx, w, nw, conv_r, conv_i);
            if (nx > 0) {
                memcpy(mix_r, conv_r + nw / 2, (size_t)(4 * nx));
                memcpy(mix_i, conv_i + nw / 2, (size_t)(4 * nx));
            }
            free(conv_r);
            free(conv_i);
            free(osc2);
            for (i = 0; i < ndec; ++i) {
                dec_r[i] = mix_r[i * decim];
                dec_i[i] = mix_i[i * decim];
            }
            x2as_decimated(dec_r, dec_i, nwork, wf, ws, nh);
            {
                float *interp = cig_interp(axis, dec_r, ndec, query, nfrm);
                for (t = 0; t < nfrm; ++t)
                    SNR2[t][b] = interp[t];
                free(interp);
            }
        }

        free(w);
        free(wf);
        free(ws);
        free(dec_r);
        free(dec_i);
        free(mix_r);
        free(mix_i);
        free(axis);

        ifdetector *ifd = cig_create_ifdetector(fc[b], fres[b]);
        ifdetector *ifd2 = cig_create_ifdetector(fc[b] + fc[b], fres[b]);
        for (t = 0; t < nfrm; ++t) {
            int nhif = ifd->nh;
            float *frame = fetch_frame(x, nx, naxis[t], nhif);
            double value;
            value = cig_ifdetector_estimate(ifd, frame, nhif);
            value = fc[b] * 0.5 < value ? value : fc[b] * 0.5;
            IF1[t][b] = fc[b] * 1.5 >= value ? value : fc[b] * 1.5;
            if (IF2) {
                value = cig_ifdetector_estimate(ifd2, frame, ifd2->nh);
                value = fc[b] * 1.5 < value ? value : fc[b] * 1.5;
                IF2[t][b] = fc[b] * 2.5 >= value ? value : fc[b] * 2.5;
            }
            free(frame);
        }
        cig_delete_ifdetector(ifd);
        cig_delete_ifdetector(ifd2);
    }

    free(query);
}
