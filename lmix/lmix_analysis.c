#include "lmix.h"
#include "lmix_utils.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <wchar.h>

#include "../config.h"
#include "../external/ciglet/ciglet.h"
#include "../external/libgvps/gvps.h"
#include "../external/libllsm/envelope.h"
#include "../external/mrq/mrq.h"
#include "../logging.h"
#include "../nebula/nebula.h"
#include "../utils.h"

// Viterbi neighborhood from context[1]
static int lmix_unwrap_fntran(void *context, int t)
{
    (void)t;
    return ((int *)context)[1];
}

// Gaussian transition weight exp(-(ds/n)^2 * 40)
static FP_TYPE lmix_unwrap_ftran(void *context, int ds, int t)
{
    (void)t;
    FP_TYPE x = (FP_TYPE)ds / (FP_TYPE)((int *)context)[0];
    return exp_2(x * -x * 40.0f);
}

// circular Viterbi unwrap of a voiced-run phase track
void lmix_unwrap_phase_track(
    FP_TYPE *phase, FP_TYPE *weight, int n, int nbin, int ntran)
{
    int nstate = 2 * nbin;
    gvps_obsrv *obs = gvps_obsrv_create_full(n, nstate);
    for (int t = 0; t < n; ++t) {
        int center = (int)((double)nbin * ((double)phase[t] / M_PI + 1.0));
        int span = sign_change_span(weight, n, t);
        for (int j = 0; j < nstate; ++j) {
            if (weight[t] <= 0.0f || span <= 9) {
                obs->slice[t]->pair[j].p =
                    0.5f / (FP_TYPE)nbin;
            } else {
                int dist = abs(center - j);
                if (nbin < dist)
                    dist = nstate - dist;
                obs->slice[t]->pair[j].p =
                    1.0f - (FP_TYPE)dist / (FP_TYPE)nbin;
                if (t == 0 || t == n - 1)
                    obs->slice[t]->pair[j].p =
                        dist != 0 ? 0.0f : 1000.0f;
            }
        }
    }

    int ctx[2];
    ctx[0] = nbin;
    ctx[1] = ntran;
    int *dst = (int *)calloc((size_t)n, sizeof(int));
    gvps_full_circular(
        dst, ctx, obs, lmix_unwrap_ftran, lmix_unwrap_fntran);
    for (int k = 0; k < n; ++k) {
        if (!(weight[k] > 0.0f &&
              sign_change_span(weight, n, k) > 9))
            continue;
        FP_TYPE target = ((FP_TYPE)dst[k] / (FP_TYPE)nbin - 1.0f) *
            (FP_TYPE)M_PI;
        FP_TYPE p = phase[k];
        if (target > p) {
            int n2 = (int)((target - p) * 0.5 / M_PI + 1.0);
            p += (FP_TYPE)M_PI * ((FP_TYPE)n2 + (FP_TYPE)n2);
        }
        FP_TYPE rem = (FP_TYPE)fmod((double)(p - target) + M_PI, 2.0 * M_PI);
        rem -= (FP_TYPE)M_PI;
        if ((FP_TYPE)M_PI / (FP_TYPE)nbin +
                (FP_TYPE)M_PI / (FP_TYPE)nbin > fabsf(rem))
            target = phase[k];
        phase[k] = target;
    }
    gvps_obsrv_free(obs);
    free(dst);
}

// mark harmonic/noise collisions and replace flagged amplitudes from a delayed moving average
static void lmix_noise_distortion_reduce_stage(llsm_layer0 *layer0)
{
    int nfrm = layer0->conf.nfrm;
    int nnos = layer0->conf.nnos;
    int nhar_max = layer0->conf.nhar;
    FP_TYPE *freqwarp = llsm_warp_freq(
        0.0f, layer0->conf.nosf, nnos, layer0->conf.noswarp);
    FP_TYPE *bin_idx = (FP_TYPE *)calloc((size_t)nnos, sizeof(FP_TYPE));
    FP_TYPE *har_freq = (FP_TYPE *)calloc((size_t)nhar_max, sizeof(FP_TYPE));
    unsigned char **flags =
        (unsigned char **)calloc((size_t)nfrm, sizeof(unsigned char *));
    if (nfrm != 0) {
        for (int i = 0; i < nfrm; ++i)
            flags[i] = (unsigned char *)calloc(12, 1);
    }
    if (nnos > 0) {
        for (int k = 0; k < nnos; ++k)
            bin_idx[k] = (FP_TYPE)k;
    }

    FP_TYPE *f0_track = (FP_TYPE *)calloc((size_t)nfrm, sizeof(FP_TYPE));
    for (int i = 0; i < nfrm; ++i) {
        llsm_frame *frame = layer0->frames[i];
        f0_track[i] = frame->f0;
        int nhar = frame->sinu->nhar;
        if (nhar <= 0)
            continue;
        int nuse = nhar <= 12 ? nhar : 12;
        for (int k = 0; k < nuse; ++k)
            har_freq[k] = ((FP_TYPE)k + 1.0f) * frame->f0;
        FP_TYPE *interp = cig_interp(freqwarp, bin_idx, nnos, har_freq, nuse);
        for (int k = 0; k < nuse; ++k) {
            int idx = (int)round(interp[k]);
            int nsum = 0;
            FP_TYPE sum = 0.0f;
            for (int d = -2; d <= 2; ++d) {
                int b = idx + d;

                if (b >= 0 && b < nnos) {
                    ++nsum;
                    sum += exp_2(frame->noise->spec[b]);
                }
            }
            if (nsum > 0 &&
                sum / (FP_TYPE)nsum > 0.3f * frame->sinu->ampl[k]) {
                if (i - 2 >= 0)
                    flags[i - 2][k] = 1;
                if (i - 1 >= 0)
                    flags[i - 1][k] = 1;
                flags[i][k] = 1;
                if (i + 1 < nfrm)
                    flags[i + 1][k] = 1;
                if (i + 2 < nfrm)
                    flags[i + 2][k] = 1;
            }
        }
        free(interp);
    }
    free(freqwarp);
    free(bin_idx);
    free(har_freq);

    FP_TYPE *series = (FP_TYPE *)calloc((size_t)nfrm, sizeof(FP_TYPE));
    for (int harm = 0; harm != 12; ++harm) {
        FP_TYPE *avg;
        if (nfrm <= 0) {
            avg = cig_moving_avg(series, nfrm, 5.0f);
        } else {
            for (int i = 0; i < nfrm; ++i) {
                int nuse = layer0->frames[i]->sinu->nhar;
                if (nuse > 12)
                    nuse = 12;
                if (harm >= nuse)
                    series[i] = 0.0f;
                else
                    series[i] =
                        layer0->frames[i]->sinu->ampl[harm];
            }
            avg = cig_moving_avg(series, nfrm, 5.0f);
            for (int i = 0; i < nfrm; ++i) {
                if (flags[i][harm] != 1)
                    continue;
                llsm_sinframe *sinu = layer0->frames[i]->sinu;
                if (harm >= sinu->nhar)
                    continue;
                int dist = sign_change_span(f0_track, nfrm, i);
                if (dist <= 4)
                    continue;
                sinu->ampl[harm] = avg[i + 2]; // no bounds check
            }
        }
        free(avg);
    }
    free(series);
    if (nfrm != 0) {
        for (int i = 0; i < nfrm; ++i)
            free(flags[i]);
    }
    free(flags);
    free(f0_track);
}

// unwrap harmonic phases over each voiced run
static void lmix_unwrap_harmonic_phase(
    llsm_layer1 *layer1, llsm_layer0 *layer0)
{
    int nfrm = layer0->conf.nfrm;
    int i = 0;
    if (nfrm <= 0)
        return;
    do {
        int j = i;
        if (layer0->frames[i]->f0 > 0.0f) {
            do {
                if (nfrm <= ++j)
                    break;
            } while (layer0->frames[j]->f0 > 0.0f);

            int n = j - i;
            FP_TYPE *phase = (FP_TYPE *)calloc((size_t)n, sizeof(FP_TYPE));
            FP_TYPE *weight = (FP_TYPE *)calloc((size_t)n, sizeof(FP_TYPE));
            int nhar_est = (int)round(7000.0f / layer0->frames[i]->f0);
            if (nhar_est > 0) {
                for (int h = 0; h < nhar_est; ++h) {
                    if (j <= i) {
                        lmix_unwrap_phase_track(
                            phase, weight, n, 30, 2);
                        continue;
                    }
                    for (int t = i; t < j; ++t) {
                        int off = t - i;
                        if (layer0->frames[t]->sinu->nhar <= h) {
                            phase[off] = 0.0f;
                            weight[off] = 0.0f;
                        } else {
                            FP_TYPE p =
                                layer1->vs_har_phse[t][h];
                            int nwrap = (int)round(
                                p * 0.1591549430918953);
                            weight[off] = 1.0f;
                            phase[off] = p -
                                (FP_TYPE)nwrap *
                                6.283185307179586f;
                        }
                    }
                    lmix_unwrap_phase_track(
                        phase, weight, n, 30, 2);
                    for (int t = i; t < j; ++t) {
                        if (layer0->frames[t]->sinu->nhar > h)
                            layer1->vs_har_phse[t][h] =
                                phase[t - i];
                    }
                }
            }
            free(phase);
            free(weight);
        }
        i = j + 1;
    } while (nfrm > i);
}

// subtract mean tail-frame noise energy from every frame
static void lmix_normalize_noise_energy(llsm_layer0 *layer0)
{
    int nfrm = layer0->conf.nfrm;
    int last_voiced = 0;
    if (nfrm > 0) {
        for (int i = 0; i < nfrm; ++i) {
            if (layer0->frames[i]->f0 > 0.0f)
                last_voiced = i;
        }
    }
    int mid = (nfrm + last_voiced) / 2;
    int n_tail = nfrm - mid;
    if (n_tail <= 2)
        return;

    int nnos = layer0->conf.nnos;
    FP_TYPE *acc = (FP_TYPE *)calloc((size_t)nnos, sizeof(FP_TYPE));
    if (nfrm > mid && nnos > 0) {
        for (int i = mid; i < nfrm; ++i) {
            FP_TYPE *spec = layer0->frames[i]->noise->spec;
            for (int j = 0; j < nnos; ++j)
                acc[j] += exp_2(spec[j]);
        }
        for (int j = 0; j < nnos; ++j)
            acc[j] *= 1.0f / (FP_TYPE)n_tail;
        for (int i = 0; i < nfrm; ++i) {
            FP_TYPE *spec = layer0->frames[i]->noise->spec;
            for (int j = 0; j < nnos; ++j) {
                FP_TYPE v = exp_2(spec[j]) - acc[j];
                if (v < 1.0e-10f)
                    v = 1.0e-10f;
                spec[j] = log_2(v);
            }
        }
    }
    free(acc);
}

// PCM + F0 contour to LMIX via layer0/layer1 analysis
lmix_data *lmix_from_wave(
    FP_TYPE *samples, int nsample, int fs,
    FP_TYPE *f0, int nf0, int nhop)
{
    log_file_printf("lmix_from_wave: initializing analysis...\n");
    llsm_parameters param = llsm_init(4);
    param.a_nnos = 96;
    param.a_tfft = 0.04f;
    param.a_nosbandf[0] = g_analysis_llsm_mbn_band0;
    param.a_nosbandf[1] = g_analysis_llsm_mbn_band1;
    param.a_nosbandf[2] = g_analysis_llsm_mbn_band2;
    param.a_nhop = nhop;
    param.a_mvf = g_analysis_llsm_max_voiced_frequency;
    param.a_nhar = g_analysis_llsm_max_num_harmonics;
    param.a_nosf = g_analysis_llsm_max_noise_frequency;
    param.a_nhare = g_analysis_llsm_energy_envelope_num_harmonics;

    log_file_printf("lmix_from_wave: lv0 LLSM analysis...\n");
    llsm_layer0 *layer0 =
        llsm_layer0_analyze(param, samples, nsample, fs, f0, nf0, NULL);

    log_file_printf("lmix_from_wave: lv1 LLSM analysis...\n");
    llsm_layer1 *layer1 =
        llsm_layer1_from_layer0(param, layer0, g_lmix_nfft, fs);

    if (g_analysis_anti_distortion != 0) {
        log_file_printf(
            "lmix_from_wave: reducing noise distortion (stage 1)...\n");
        lmix_noise_distortion_reduce_stage(layer0);
        log_file_printf(
            "lmix_from_wave: reducing noise distortion (stage 2)...\n");
        lmix_unwrap_harmonic_phase(layer1, layer0);
    }
    if (g_analysis_noise_reduction != 0) {
        log_file_printf("lmix_from_wave: reducing noise...\n");
        lmix_normalize_noise_energy(layer0);
    }

    log_file_printf(
        "lmix_from_wave: restructuring data (%d frames)...\n", nf0);
    FP_TYPE *faxis = llsm_uniform_faxis(g_lmix_nfft, fs);
    int ns = g_lmix_nfft / 2 + 1;
    FP_TYPE faxis_upper = faxis[ns - 1] + faxis[ns - 1] - faxis[ns - 2];
    FP_TYPE *har_freq = (FP_TYPE *)calloc(
        (size_t)g_lmix_nband_default * (size_t)layer0->conf.nhar,
        sizeof(FP_TYPE));
    lmix_data *lmix = lmix_create(layer0->conf);

    int i = 0;
    if (nf0 > 0) {
        do {
            llsm_frame *lv0 = layer0->frames[i];
            lmix->frames[i] = lmix_create_frame(
                layer0->conf, g_lmix_nfft, lv0->sinu->nhar);
            lmix_frame *frame = lmix->frames[i];
            int nband = g_lmix_nband_default;
            frame->vt_stride = nband;
            int nhar = frame->nhar;
            frame->f0 = lv0->f0;
            frame->rd = layer1->vs_rd[i];
            frame->vt_freq_step = lv0->f0 / (FP_TYPE)nband;

            if (lv0->wave != NULL) {
                if (frame->wave == NULL) {
                    frame->wave = llsm_create_wave(
                        lv0->wave->samples,
                        lv0->wave->nsample,
                        lv0->wave->fs);
                }
                llsm_copy_wave(frame->wave, lv0->wave);
                frame = lmix->frames[i];
                lv0 = layer0->frames[i];
                nhar = frame->nhar;
            }

            llsm_nosframe *src_noise = lv0->noise;
            llsm_nosframe *dst_noise = frame->noise;
            dst_noise->nchannel = src_noise->nchannel;
            dst_noise->nnos = src_noise->nnos;
            if (src_noise->nnos > 0) {
                for (int k = 0; k != src_noise->nnos; ++k)
                    dst_noise->spec[k] = src_noise->spec[k];
            }
            if (src_noise->nchannel > 0) {
                for (int ch = 0; ch != src_noise->nchannel; ++ch) {
                    dst_noise->emin[ch] = src_noise->emin[ch];
                    dst_noise->eenv[ch]->nhar = 0;
                }
            }

            if (nhar > 0) {
                for (int k = 0; k != nhar; ++k)
                    frame->vs_har_phse[k] = layer1->vs_har_phse[i][k];
                frame->vs_har_phse[0] = 0.0f;
                int nspec = nhar * nband;
                free(frame->vt_magn);
                if (nspec > 0) {
                    for (int k = 0; k != nspec; ++k)
                        har_freq[k] = ((FP_TYPE)k + 1.0f) * frame->vt_freq_step;
                }
                frame->vt_magn = cig_interpu(
                    faxis[0], faxis_upper, layer1->vt_resp_magn[i], ns,
                    har_freq, nspec);

                if (src_noise->nchannel > 0) {
                    for (int ch = 0; ch != src_noise->nchannel; ++ch) {
                        llsm_sinframe *dst_eenv = dst_noise->eenv[ch];
                        llsm_sinframe *src_eenv = src_noise->eenv[ch];

                        dst_eenv->nhar = src_eenv->nhar;
                        if (src_eenv->nhar > 0) {
                            FP_TYPE ph0 = lv0->sinu->phse[0];
                            for (int m = 0; m != src_eenv->nhar; ++m) {
                                dst_eenv->ampl[m] = src_eenv->ampl[m];
                                dst_eenv->phse[m] = src_eenv->phse[m] -
                                    ((FP_TYPE)m + 1.0f) * ph0;
                            }
                        }
                    }
                }
            }
            ++i;
        } while (nf0 != i);
    }

    free(faxis);
    free(har_freq);
    llsm_delete_layer0(layer0);
    llsm_delete_layer1(layer1);
    llsm_deinit(param);
    log_file_printf("lmix_from_wave: analysis finished.\n");
    return lmix;
}

// read every second float from an UTAU _wav.frq file
static FP_TYPE *lmix_load_frq_file(
    const wchar_t *path, FP_TYPE *hop_out, int *count_out)
{
    FILE *stream = _wfopen(path, L"rb");
    if (stream == NULL)
        return NULL;

    puts("load_frq_file: frq file detected. Loading...");
    unsigned char scratch[36];
    int hop = 0;
    if (fread(scratch, 8, 1, stream) != 1 ||
        fread(&hop, 4, 1, stream) != 1 ||
        fread(scratch, 4, 1, stream) != 1) {
        fclose(stream);
        return NULL;
    }

    long data_begin = ftell(stream);
    if (data_begin < 0 || fseek(stream, 0, SEEK_END) != 0) {
        fclose(stream);
        return NULL;
    }
    long data_end = ftell(stream);
    if (data_end < data_begin || fseek(stream, data_begin, SEEK_SET) != 0) {
        fclose(stream);
        return NULL;
    }

    unsigned int count = (unsigned int)(data_end - data_begin) >> 3;
    FP_TYPE *f0 = (FP_TYPE *)calloc((size_t)count, sizeof(*f0));
    if (count != 0 && f0 == NULL) {
        fclose(stream);
        return NULL;
    }

    for (unsigned int i = 0; i < count; ++i) {
        if (fread(scratch, 4, 1, stream) != 1 ||
            fread(&f0[i], 4, 1, stream) != 1) {
            free(f0);
            fclose(stream);
            return NULL;
        }
    }
    fclose(stream);

    *hop_out = (FP_TYPE)hop;
    *count_out = (int)count;
    return f0;
}

// load a wav as a single-wave lmix without spectral analysis
lmix_data *lmix_from_file_direct(const wchar_t *path)
{
    int wav_fs = 0;
    int wav_bits = 0;
    int wav_nsample = 0;
    log_file_printf("lmix_from_file_direct: loading audio...\n");
    FP_TYPE *wav = audioread(path, &wav_fs, &wav_bits, &wav_nsample);
    if (wav == NULL) {
        log_tee_printf("Error: audioread failed for %ls.\n", path);
        return NULL;
    }

    int hop = hop_from_fs(g_output_sampling_rate);
    FP_TYPE thop = (FP_TYPE)hop / (FP_TYPE)g_output_sampling_rate;
    double frame_count = (double)wav_nsample / ((double)wav_fs * thop);
    int nframe = (int)round(frame_count);

    FP_TYPE *nosbandf = (FP_TYPE *)calloc(3, sizeof(FP_TYPE));
    nosbandf[0] = g_analysis_llsm_mbn_band0;
    nosbandf[1] = g_analysis_llsm_mbn_band1;
    nosbandf[2] = g_analysis_llsm_mbn_band2;

    llsm_conf conf;
    memset(&conf, 0, sizeof(conf));
    conf.nfrm = nframe;
    conf.nhop = hop;
    conf.nhar = g_analysis_llsm_max_num_harmonics;
    conf.nhare = g_analysis_llsm_energy_envelope_num_harmonics;
    conf.nnos = 96;
    conf.noswarp = 15000.0f;
    conf.mvf = g_analysis_llsm_max_voiced_frequency;
    conf.nosf = g_analysis_llsm_max_noise_frequency;
    conf.thop = (FP_TYPE)thop;
    conf.nosbandf = nosbandf;
    conf.nnosband = 4;

    lmix_data *lmix = lmix_create(conf);
    for (int i = 0; i < nframe; ++i)
        lmix->frames[i] = lmix_create_frame(conf, g_lmix_nfft, 0);

    lmix->frames[0]->wave = llsm_create_wave(wav, wav_nsample, wav_fs);
    lmix->frames[0]->wave->t0 = 0.0f;
    free(nosbandf);
    free(wav);
    return lmix;
}

// load from .llsm / desc.mrq cache, or analyze the wav; skip_llsm is flag _f
lmix_data *lmix_from_file(
    const wchar_t *path, float start, float end, int skip_llsm)
{
    size_t path_len = wcslen(path);
    wchar_t *llsm_path = wcscat(wcscpy(
        (wchar_t *)calloc(path_len + 6, sizeof(wchar_t)), path), L".llsm");
    size_t dir_len = path_dir_prefix_len(path);
    wchar_t *mrq_path = wcscat(wcsncpy(
        (wchar_t *)calloc(path_len + 8, sizeof(wchar_t)), path, dir_len),
        L"desc.mrq");

    time_t wav_mtime = file_mtime(path);
    xchg_node *root = xchg_import(llsm_path, start, end);
    lmix_data *lmix = NULL;

    if (root == NULL) {
        log_tee_printf(
            "Data record is not detected. The required file will be generated in the current session.\n");
    } else if (g_auto_update_llsm_mrq &&
               wav_mtime > file_mtime(llsm_path)) {
        log_tee_printf(
            "The .wav file is newer than the data record. The record will be updated in the current session.\n");
    } else {
        xchg_node *data = xchg_node_get_field(root, "data");
        xchg_node *version_node = xchg_node_get_field(data, "version");
        if (version_node == NULL) {
            log_tee_printf(
                "An invalid data record is detected. This file will be re-generated in the current session.\n");
        } else {
            int version = (int)version_node->float_value;
            if ((unsigned int)(version - 520) > 0x4FU) {
                if (version > 519) {
                    log_tee_printf(
                        "A data record with newer but incompatible version is detected. Synthesis will terminate for overwrite-protection. You may consider either update moresampler or delete the data record which will be automatically re-generated by next session.\n");
                    return NULL;
                }
                log_tee_printf(
                    "A data record with older and incompatible version is detected. This record will be updated in the current session.\n");
                if (data == NULL) {
                    log_tee_printf(
                        "An invalid data record is detected. This file will be re-generated in the current session.\n");
                }
            } else {
                if (version > 524)
                    log_tee_printf(
                        "A data record with newer version is detected. You may consider update moresampler.\n");
                else if (version != 524)
                    log_tee_printf(
                        "A data record with older version is detected. You may consider delete the data record and moresampler will automatically re-analyze by next session.\n");
                if (data == NULL) {
                    log_tee_printf(
                        "An invalid data record is detected. This file will be re-generated in the current session.\n");
                } else {
                    lmix = lmix_from_xchg(xchg_node_get_field(root, "data"));
                    if (lmix == NULL) {
                        log_tee_printf(
                            "An invalid data record is detected. This file will be re-generated in the current session.\n");
                        xchg_node_destroy(root);
                    } else {
                        xchg_node_destroy(root);
                        free(mrq_path);
                        free(llsm_path);
                        return lmix;
                    }
                }
            }
        }
    }

    int wav_fs = 0;
    int wav_bits = 0;
    int wav_nsample = 0;
    log_file_printf("lmix_from_file: loading audio...\n");
    FP_TYPE *wav = audioread(path, &wav_fs, &wav_bits, &wav_nsample);
    if (wav == NULL) {
        log_tee_printf("Error: audioread failed for %ls.\n", path);
        return NULL;
    }

    mrq_entry *entry = NULL;
    FP_TYPE *frq = NULL;
    FP_TYPE frq_hop = 0.0f;
    int frq_count = 0;
    int hop = hop_from_fs(wav_fs);
    int f0_hop = hop;
    FP_TYPE *f0 = NULL;
    int nf0 = 0;

    log_file_printf("lmix_from_file: loading mrq...\n");
    FILE *mrq = mrq_open(mrq_path, L"rw");
    if (mrq == NULL) {
        log_tee_printf("Warning: cannot load or create mrq file.\n");
    } else {
        log_file_printf("lmix_from_file: locking mrq...\n");
        xchg_file_lock(mrq);
        log_file_printf("lmix_from_file: getting mrq version...\n");
        if ((unsigned int)(mrq_get_version(mrq) - 1) > 8U) {
            fclose(mrq);
            log_tee_printf(
                "A mrq data descriptor of unsupported version is detected. The file is ignored.\n");
        }
        // still calls mrq_get_entry / fclose after the unsupported-version branch
        entry = mrq_get_entry(mrq, path_get_filename((wchar_t *)path));
        fclose(mrq);
    }

    if (entry != NULL && g_auto_update_llsm_mrq) {
        if (entry->timestamp != 0 && wav_mtime > (time_t)entry->timestamp) {
            log_tee_printf(
                "A mrq data entry is detected but it is older than the .wav file. Moresampler will update the entry in current session.\n");
            entry = NULL;
        } else if (entry->timestamp != 0 &&
                   entry->modified != 0 && entry->modified <= 519) {
            log_tee_printf(
                "A mrq data entry is detected but it is generated by a deprecated version of Moresampler. Moresampler will update the entry in current session.\n");
            entry = NULL;
        }
    }

    if (entry != NULL) {
        log_tee_printf("A mrq data entry is found and loaded.\n");
        nf0 = entry->nf0;
        f0_hop = entry->nhop;
        f0 = (FP_TYPE *)calloc((size_t)nf0, sizeof(*f0));
        for (int i = 0; i < nf0; ++i)
            f0[i] = (FP_TYPE)entry->f0[i];
        hop = hop_from_fs(wav_fs);
    } else {
        if (g_load_frq > 0) {
            wchar_t *frq_path = wcscpy(
                (wchar_t *)calloc(path_len + 5, sizeof(wchar_t)), path);
            frq_path[wcslen(frq_path) - 4] = L'_';
            frq_path = wcscat(frq_path, L".frq");
            frq = lmix_load_frq_file(frq_path, &frq_hop, &frq_count);
            if (frq != NULL)
                log_file_printf("A frq file is found and loaded.\n");
            free(frq_path);
        }

        FP_TYPE path_f0 = f0_from_dirname(path);
        hop = hop_from_fs(wav_fs);
        log_file_printf("f0_from_wave: initializing analysis...\n");
        FP_TYPE fmin = g_analysis_f0_min;
        FP_TYPE fmax = g_analysis_f0_max;
        if (g_analysis_f0_range_from_path && path_f0 != 0.0f) {
            log_tee_printf(
                "f0_from_wave: average f0 estimated from directory name = %fHz\n",
                path_f0);
            fmin = 0.6f * path_f0;
            fmax = 1.7f * path_f0;
        }
        int biased = g_analysis_biased_f0_estimation != 0;
        FP_TYPE dither = biased ? 0.05f : 0.1f;
        if ((unsigned int)g_load_frq <= 1U) {
            log_file_printf("f0_from_wave: estimating pitch...\n");
            f0 = f0_from_wave(
                fmin, fmax, (FP_TYPE)hop / (FP_TYPE)wav_fs,
                1, biased, dither, 0.2f, 2.0f, 1.5f, 0.0f, 0.0f,
                3, 10, 10, wav, wav_nsample, wav_fs, &nf0);
        }
        if (frq != NULL) {
            if (g_load_frq == 1 && f0 != NULL) {
                for (int i = 0; i < nf0; ++i) {
                    if (f0[i] <= 0.0f) {
                        FP_TYPE value = frq[frq_index(
                            i * hop, frq_hop, frq_count)];
                        if (value > 0.0f)
                            f0[i] = value;
                    }
                }
            } else if (g_load_frq == 2) {
                nf0 = wav_nsample / hop;
                f0 = (FP_TYPE *)calloc((size_t)nf0, sizeof(*f0));
                for (int i = 0; i < nf0; ++i)
                    f0[i] = frq[frq_index(i * hop, frq_hop, frq_count)];
            }
        }
        entry = NULL;
        f0_hop = hop;
    }

    if (hop != f0_hop) {
        int new_count = nf0 * f0_hop / hop;
        FP_TYPE *resampled =
            (FP_TYPE *)calloc((size_t)new_count, sizeof(*resampled));

        for (int i = 0; i < new_count; ++i)
            resampled[i] = f0[i * hop / f0_hop];
        free(f0);
        f0 = resampled;
        nf0 = new_count;
        f0_hop = hop;
    }

    lmix = NULL;
    if (skip_llsm == 0) {
        lmix = lmix_from_wave(wav, wav_nsample, wav_fs, f0, nf0, f0_hop);
        lmix->duration = (FP_TYPE)wav_nsample / (FP_TYPE)wav_fs;
    }

    if (entry == NULL) {
        entry = create_mrq_entry(nf0);
        entry->nhop = f0_hop;
        entry->fs = wav_fs;
        entry->modified = 524;
        for (int i = 0; i < nf0; ++i)
            entry->f0[i] = f0[i];
        mrq = mrq_open(mrq_path, L"rw");
        if (mrq != NULL) {
            xchg_file_lock(mrq);
            if (mrq_write_entry(
                    mrq, path_get_filename((wchar_t *)path), entry) != 0)
                log_file_printf("A new mrq data entry is created.\n");
            else
                log_tee_printf("Warning: cannot write back to mrq file.\n");
            fclose(mrq);
        } else {
            log_tee_printf("Warning: cannot write back to mrq file.\n");
        }
    }

    free(f0);
    free(frq);
    delete_mrq_entry(entry);
    if (lmix != NULL) {
        log_file_printf("lmix_from_file: exporting to llsm file...\n");
        xchg_node *node = lmix_to_xchg(lmix);
        xchg_export_counted(node, llsm_path);
        xchg_node_destroy(node);
    }
    free(wav);
    free(mrq_path);
    free(llsm_path);
    return lmix;
}
