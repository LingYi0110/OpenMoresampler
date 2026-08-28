#include "nebula.h"

#include "../external/fastapprox/fasttrig.h"
#include "../external/ciglet/ciglet.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


// IF/SNR: heterodyne each band, two centred complex FIRs, residual power, interpolate
static void get_wcoef(const char *name, float a[4])
{
    a[0] = a[1] = a[2] = a[3] = 0.0f;
    if (!name)
        return;
    if (strcmp(name, "nuttall83") == 0) {
        a[0] = 0.338946f; a[1] = 0.481973f; a[2] = 0.161054f; a[3] = 0.018027f;
    } else if (strcmp(name, "nuttall98") == 0) {
        a[0] = 0.3635819f; a[1] = 0.4891775f; a[2] = 0.1365995f; a[3] = 0.0106411f;
    } else if (strcmp(name, "nuttall93") == 0) {
        a[0] = 0.355768f; a[1] = 0.487396f; a[2] = 0.144232f; a[3] = 0.012604f;
    } else if (strcmp(name, "nuttall64") == 0) {
        a[0] = 0.40897f; a[1] = 0.5f; a[2] = 0.09103f;
    } else if (strcmp(name, "hanning") == 0) {
        a[0] = 0.5f; a[1] = 0.5f;
    } else if (strcmp(name, "hamming") == 0) {
        a[0] = 0.54f; a[1] = 0.46f;
    } else if (strcmp(name, "blackman") == 0) {
        a[0] = 0.42f; a[1] = 0.5f; a[2] = 0.08f;
    } else if (strcmp(name, "blackman-harris") == 0) {
        a[0] = 0.4243801f; a[1] = 0.4973406f; a[2] = 0.0782793f;
    } else if (strcmp(name, "boxcar") == 0) {
        a[0] = 1.0f;
    }
}

// cosine-sum window, normalized to unit sum
static float *make_window(const float aw[4], int nw, float omegaw)
{
    float *w = (float *)calloc((size_t)nw, sizeof(*w));
    int half = nw / 2;

    for (int k = 0; k < 4 && aw[k] > 0.0f; ++k)
        for (int t = -half; t < nw - half; ++t)
            w[t + half] += aw[k] * fastcosfull(t * (k * omegaw));

    float sum = 0.0f;
    for (int i = 0; i < nw; ++i)
        sum += w[i];
    sum = 1.0f / sum;
    for (int i = 0; i < nw; ++i)
        w[i] *= sum;
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

static int fft_block_size(int nh)
{
    if (nh <= 0)
        return 0;
    int exp2 = (int)rintf(log2f((float)nh)) + 2;
    if (exp2 < 2)
        exp2 = 2;
    if (exp2 > 29)
        return 0;
    return 1 << exp2;
}

// complex overlap-add convolution; output length nx + nh - 1
static int fftconv_complex(const float *xr, const float *xi, int nx,
                           const float *hr, const float *hi, int nh,
                           float *yr, float *yi)
{
    int nfft = fft_block_size(nh);
    if (!xr || !hr || !yr || nx <= 0 || nh <= 0 || nfft <= nh)
        return 0;
    int stride = nfft - nh;
    int nblock = nx / stride + 1;
    float *fr = (float *)calloc((size_t)nfft, sizeof(*fr));
    float *fi = (float *)calloc((size_t)nfft, sizeof(*fi));
    float *br = (float *)calloc((size_t)nfft, sizeof(*br));
    float *bi = (float *)calloc((size_t)nfft, sizeof(*bi));
    float *work = (float *)malloc((size_t)(2 * nfft) * sizeof(*work));
    if (!fr || !fi || !br || !bi || !work) {
        free(fr); free(fi); free(br); free(bi); free(work);
        return 0;
    }

    memcpy(fr, hr, (size_t)nh * sizeof(*fr));
    if (hi)
        memcpy(fi, hi, (size_t)nh * sizeof(*fi));
    cig_fft(fr, fi, fr, fi, nfft, work, -1.0f);

    memset(yr, 0, (size_t)(nx + nh - 1) * sizeof(*yr));
    if (yi)
        memset(yi, 0, (size_t)(nx + nh - 1) * sizeof(*yi));

    for (int block = 0, offset = 0; block < nblock;
         ++block, offset += stride) {
        int copy = nx - offset;
        if (copy < 0)
            copy = 0;
        if (copy > stride)
            copy = stride;
        memset(br, 0, (size_t)nfft * sizeof(*br));
        memset(bi, 0, (size_t)nfft * sizeof(*bi));
        if (copy > 0) {
            memcpy(br, xr + offset, (size_t)copy * sizeof(*br));
            if (xi)
                memcpy(bi, xi + offset, (size_t)copy * sizeof(*bi));
        }
        cig_fft(br, bi, br, bi, nfft, work, -1.0f);
        for (int i = 0; i < nfft; ++i) {
            float ar = br[i];
            float ai = bi[i];
            br[i] = ar * fr[i] - ai * fi[i];
            bi[i] = ar * fi[i] + ai * fr[i];
        }
        cig_fft(br, bi, br, bi, nfft, work, 1.0f);
        for (int i = 0; i < nfft && offset + i < nx + nh - 1; ++i) {
            yr[offset + i] += br[i];
            if (yi)
                yi[offset + i] += bi[i];
        }
    }

    free(fr); free(fi); free(br); free(bi); free(work);
    return 1;
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
    if (!fftconv_complex(yr, yi, nwork, w, NULL, nh, y2r, y2i)) {
        free(w); free(ws); free(y2r); free(y2i); free(power); free(smooth);
        return 0;
    }
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
    if (!fftconv_complex(power, NULL, nwork, ws, NULL, nh, smooth, NULL)) {
        free(w); free(ws); free(y2r); free(y2i); free(power); free(smooth);
        return 0;
    }

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
    if (!fftconv_complex(mix_r, mix_i, nx, w, NULL, nw, conv_r, conv_i)) {
        free(w); free(osc); free(mix_r); free(mix_i);
        free(conv_r); free(conv_i); free(dec_r); free(dec_i);
        return 0;
    }
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
