#include "lmix_utils.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <wchar.h>

#include <../external/libaiff/libaiff/libaiff.h>

#include "external/ciglet/ciglet.h"
#include "lmix_lf.h"
#include "../utils.h"

int sign_change_span(FP_TYPE *x, int n, int i)
{
    int hi = i + 10;
    if (hi > n - 1)
        hi = n - 1;
    int fwd = 10;
    for (int j = i; j < hi; ++j) {
        if ((x[j] > 0.0f && x[j + 1] <= 0.0f) ||
            (x[j] <= 0.0f && x[j + 1] > 0.0f)) {
            fwd = j - i;
            break;
        }
    }

    int lo = i - 10;
    if (lo <= 0)
        lo = 1;
    int bwd = 10;
    for (int j = i; j > lo; --j) {
        if ((x[j] > 0.0f && x[j - 1] <= 0.0f) ||
            (x[j] <= 0.0f && x[j - 1] > 0.0f)) {
            bwd = i - j;
            break;
        }
    }
    return fwd <= bwd ? fwd : bwd;
}

int hop_from_fs(int fs)
{
    return (int)pow(2.0, ceil(log2((double)fs * 0.005)));
}

// Convert a note name to a MIDI note number.
int midi_from_note(const char *note)
{
    if (note == NULL)
        return 0;
    size_t n = strlen(note);
    if (n != 2 && n != 3)
        return 0;

    char letter = note[0];
    if (letter >= 'a' && letter <= 'g')
        letter = (char)(letter - ('a' - 'A'));
    if (letter < 'A' || letter > 'G')
        return 0;
    if (note[n - 1] < '0' || note[n - 1] > '9')
        return 0;

    int accidental = 0;
    if (n == 3) {
        if (note[1] == '#')
            accidental = 1;
        else if (note[1] == 'b')
            accidental = -1;
        else
            return 0;
    }

    int octave = note[n - 1] - '0';
    static const int pitch_class[7] = {9, 11, 0, 2, 4, 5, 7};
    return pitch_class[letter - 'A'] + 12 * (octave + 1) + accidental;
}

// Convert a MIDI note number to frequency in Hz.
FP_TYPE midi_to_hz(int midi)
{
    return (FP_TYPE)(pow(2.0, ((double)midi - 69.0) * (1.0 / 12.0)) *
                     440.0);
}

// Convert frequency in Hz to the Mel scale.
FP_TYPE hz_to_mel(FP_TYPE hz)
{
    return (FP_TYPE)(log_2(0.004545454545454545 * hz) * 1731.234049066752
                     + 5700.0);
}

int parse_note_f0(const char *note, FP_TYPE *f0_out)
{
    int midi = midi_from_note(note);
    if (midi == 0)
        return 0;
    *f0_out = midi_to_hz(midi);
    return 1;
}

// Infer an average F0 from a parent directory note suffix.
FP_TYPE f0_from_dirname(const wchar_t *path)
{
    if (path == NULL)
        return 0.0f;
    wchar_t *parent = _wcsdup(path);
    if (parent == NULL)
        return 0.0f;

    wchar_t *p = parent + wcslen(parent);
    while (p != parent && p[-1] != L'/' && p[-1] != L'\\')
        --p;
    if (p == parent) {
        free(parent);
        return 0.0f;
    }
    p[-1] = L'\0';

    wchar_t *folder = parent + wcslen(parent);
    while (folder != parent && folder[-1] != L'/' && folder[-1] != L'\\')
        --folder;
    size_t n = wcslen(folder);
    char *narrow = (char *)calloc(n + 1, 1);
    if (narrow == NULL) {
        free(parent);
        return 0.0f;
    }
    for (size_t i = 0; i < n; ++i)
        narrow[i] = (char)folder[i];

    FP_TYPE f0;
    if (parse_note_f0(narrow, &f0)) {
        free(narrow);
        free(parent);
        return f0;
    }
    for (size_t i = 0; i < n; ++i) {
        if (narrow[i] == '_' && parse_note_f0(narrow + i + 1, &f0)) {
            free(narrow);
            free(parent);
            return f0;
        }
    }

    free(narrow);
    free(parent);
    return 0.0f;
}

// Return the number of wchar_t elements through the last path separator.
size_t path_dir_prefix_len(const wchar_t *path)
{
    if (path == NULL)
        return 0;
    size_t n = wcslen(path);
    while (n > 0) {
        wchar_t c = path[n - 1];

        if (c == L'/' || c == L'\\')
            return n;
        --n;
    }
    return 0;
}

// Return st_mtime even if _wstat fails.
time_t file_mtime(const wchar_t *path)
{
#if MORESAMPLER_ARCH == 64
    struct _stat64 st;

    _wstat64(path, &st);
#else
    struct _stat32 st;

    _wstat32(path, &st);
#endif
    return (time_t)st.st_mtime;
}

// Convert a sample position to a frequency-hop index, clamped to count.
int frq_index(int sample, FP_TYPE frq_hop, int count)
{
    double index = (double)sample * (1.0 / (double)frq_hop);
    if ((double)count < index)
        index = (double)count;
    return (int)index;
}

FP_TYPE lerp_phase(FP_TYPE a, FP_TYPE b, FP_TYPE t)
{
    return (FP_TYPE)atan2(
        (double)sin_2(b) + ((double)sin_2(a) - sin_2(b)) * t,
        (double)cos_2(b) + ((double)cos_2(a) - cos_2(b)) * t);
}

// Read .wav via wavread or AIFF via libaiff; multi-channel AIFF keeps channel 0.
FP_TYPE *audioread(
    const wchar_t *path, int *fs, int *nbit, int *nsample)
{
    signed int n = (signed int)wcslen(path);
    char *mbs = wchar_to_mbs(path);
    if (n > 3) {
        const wchar_t *last4 = path + (n - 4);

        if (wcscmp(last4, L".wav") == 0) {
            FP_TYPE *wav = wavread(mbs, fs, nbit, nsample);
            free(mbs);
            return wav;
        }
        if (wcscmp(last4, L".aif") == 0
            || (n != 4
                && (wcscmp(path + (n - 5), L".aiff") == 0
                    || wcscmp(path + (n - 5), L".aifc") == 0))) {
            AIFF_Ref ref = AIFF_OpenFile(mbs, F_RDONLY);
            free(mbs);
            if (ref == NULL)
                return NULL;
            uint64_t nframes;
            int channels;
            double rate;
            int segment;
            if (AIFF_GetAudioFormat(ref, &nframes, &channels, &rate,
                    nbit, &segment) > 0) {
                *fs = (int)rate;
                *nsample = (int)nframes;
                float *interleaved = (float *)calloc(
                    (size_t)channels * (size_t)nframes,
                    sizeof(float));
                AIFF_Seek(ref, 0);
                if (AIFF_ReadSamplesFloat(ref, interleaved,
                        channels * (int)nframes) >= 0) {
                    float *mono = (float *)calloc(
                        (size_t)nframes, sizeof(float));
                    if (nframes != 0) {
                        for (size_t i = 0; i < (size_t)nframes; ++i)
                            mono[i] = interleaved[i * (size_t)channels];
                    }
                    free(interleaved);
                    AIFF_CloseFile(ref);
                    return mono;
                }
                free(interleaved);
            }
            AIFF_CloseFile(ref);
            return NULL;
        }
    }
    free(mbs);
    return NULL;
}

/* Pitch-shift frequency grid; keep its log-frequency endpoints typed. */
FP_TYPE *lmix_pitch_freqwarp(
    FP_TYPE low, FP_TYPE high, int n, FP_TYPE warp)
{
    lmix_warp_bounds bounds = lmix_log_frequency_bounds(low, high, warp);
    FP_TYPE *freq = (FP_TYPE *)calloc((size_t)n + 1, sizeof(FP_TYPE));
    for (int i = 0; i <= n; ++i)
        freq[i] = warp * (exp(((FP_TYPE)i / n *
            (bounds.high - bounds.low) + bounds.low) / 5000.0) - 1.0);
    return freq;
}

FP_TYPE *lmix_lfmodel_spectrum(
    lfmodel model, FP_TYPE *freq, int nf, FP_TYPE *dst_phase)
{
    lmix_lfparam p = lmix_lfparam_from_model(model);
    FP_TYPE e = p.e;
    FP_TYPE a = p.a;
    FP_TYPE wg = p.wg;
    FP_TYPE E0 = p.E0;
    FP_TYPE sin_wgTe = p.sin_wgTe;
    FP_TYPE cos_wgTe = p.cos_wgTe;
    FP_TYPE Te = p.Te;
    FP_TYPE Ta = p.Ta;
    FP_TYPE T0 = p.T0;
    FP_TYPE *dst_magn = (FP_TYPE *)calloc((size_t)nf, sizeof(FP_TYPE));
    FP_TYPE e1eTa = e * (1.0 - e * Ta);
    double tail_phase_scale = -2.0 * M_PI * (T0 - Te);
    double te_phase_scale = -2.0 * M_PI * Te;
    double eTa_omega_scale = 2.0 * M_PI * (e * Ta);
    for (int i = 0; i < nf; ++i) {
        double Omega = 2.0 * M_PI * freq[i];
        cplx asubipif = c_cplx(a, -Omega);
        cplx P1 = c_div(c_cplx(E0, 0.0f),
                        c_add(c_mul(asubipif, asubipif), c_cplx(wg * wg, 0.0f)));
        cplx P2 = c_add(c_cplx(wg, 0.0f),
            c_mul(c_exp_2(c_mul(asubipif, c_cplx(Te, 0.0f))),
                  c_sub(c_mul(asubipif, c_cplx(sin_wgTe, 0.0f)),
                        c_cplx(wg * cos_wgTe, 0.0f))));
        cplx P3 = c_mul(c_cplx(model.Ee, 0.0f), c_div(
            c_exp_2(c_cplx(0.0f, te_phase_scale * freq[i])),
            c_mul(c_cplx(0.0f, eTa_omega_scale * freq[i]), c_cplx(e, Omega))));
        cplx P4 = c_sub(
            c_mul(c_cplx(e1eTa, 0.0f),
                  c_sub(c_cplx(1.0f, 0.0f),
                        c_exp_2(c_cplx(0.0f, tail_phase_scale * freq[i])))),
            c_cplx(0.0f, eTa_omega_scale * freq[i]));
        cplx G = c_add(c_mul(P1, P2), c_mul(P3, P4));
        dst_magn[i] = c_abs(G);
        if (isnan(dst_magn[i])) {
            dst_magn[i] = 0.0f;
            G = c_cplx(0.0f, 0.0f);
        }
        if (dst_phase != NULL)
            dst_phase[i] = c_arg(G);
    }
    return dst_magn;
}
