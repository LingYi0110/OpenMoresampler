#include "lmix.h"

#include <stdlib.h>


int g_lmix_nfft = 2048;
int g_lmix_nband_default = 3;

lmix_data *lmix_create(llsm_conf conf)
{
    lmix_data *lmix = (lmix_data *)calloc(1, sizeof(lmix_data));
    lmix->conf = conf;

    int nfaxis = conf.nnosband - 1;
    lmix->conf.nosbandf =
        (FP_TYPE *)calloc((size_t)nfaxis, sizeof(FP_TYPE));
    if (nfaxis > 0) {
        memcpy(lmix->conf.nosbandf, conf.nosbandf,
               (size_t)nfaxis * sizeof(FP_TYPE));
    }

    lmix->version = 524;
    lmix->duration = conf.nfrm * conf.thop;
    if (conf.nfrm > 0) {
        lmix->frames = (lmix_frame **)calloc(
            (size_t)conf.nfrm, sizeof(lmix_frame *));
    } else {
        lmix->frames = NULL;
    }
    return lmix;
}

lmix_frame *lmix_create_frame(llsm_conf conf, int nfft, int nhar)
{
    (void)nfft;

    lmix_frame *frame = (lmix_frame *)calloc(1, sizeof(lmix_frame));
    frame->nhar = nhar;
    frame->f0 = 0.0f;

    if (nhar > 0) {
        int matrix_count = nhar * g_lmix_nband_default;
        frame->vt_magn = (FP_TYPE *)calloc(
            (size_t)matrix_count, sizeof(FP_TYPE));
        frame->vs_har_phse = (FP_TYPE *)calloc(
            (size_t)nhar, sizeof(FP_TYPE));

        FP_TYPE *cursor = frame->vt_magn;
        for (int i = 0; i < matrix_count; ++i)
            *cursor++ = -99.0f;
    } else {
        frame->vt_magn = NULL;
        frame->vs_har_phse = NULL;
    }

    frame->wave = NULL;

    llsm_nosframe *noise = (llsm_nosframe *)calloc(1, sizeof(llsm_nosframe));
    frame->noise = noise;
    noise->nchannel = conf.nnosband;
    noise->nnos = conf.nnos;
    noise->spec = (FP_TYPE *)calloc(
        (size_t)conf.nnos, sizeof(FP_TYPE));
    noise->emin = (FP_TYPE *)calloc(
        (size_t)conf.nnosband, sizeof(FP_TYPE));
    noise->eenv = (llsm_sinframe **)calloc(
        (size_t)conf.nnosband, sizeof(llsm_sinframe *));

    FP_TYPE *cursor = noise->spec;
    for (int i = 0; i < conf.nnos; ++i)
        *cursor++ = -100.0f;

    for (int i = 0; i < conf.nnosband; ++i) {
        llsm_sinframe *eenv =
            (llsm_sinframe *)calloc(1, sizeof(llsm_sinframe));
        noise->eenv[i] = eenv;
        eenv->ampl = (FP_TYPE *)calloc(
            (size_t)conf.nhare, sizeof(FP_TYPE));
        eenv->phse = (FP_TYPE *)calloc(
            (size_t)conf.nhare, sizeof(FP_TYPE));
        eenv->nhar = conf.nhare;
    }

    frame->growl = 0.0f;
    frame->distortion = 0.0f;
    return frame;
}

void lmix_copy_frame(lmix_frame *dst, lmix_frame *src)
{
    dst->f0 = src->f0;
    dst->vt_stride = src->vt_stride;
    dst->rd = src->rd;
    dst->vt_freq_step = src->vt_freq_step;
    dst->growl = src->growl;
    dst->distortion = src->distortion;

    if (dst->nhar < src->nhar) {
        dst->vs_har_phse = (FP_TYPE *)realloc(
            dst->vs_har_phse, (size_t)src->nhar * sizeof(FP_TYPE));
        dst->vt_magn = (FP_TYPE *)realloc(
            dst->vt_magn,
            (size_t)dst->vt_stride * (size_t)src->nhar * sizeof(FP_TYPE));
    }

    dst->nhar = src->nhar;
    for (int i = 0; i < src->nhar; ++i)
        dst->vs_har_phse[i] = src->vs_har_phse[i];

    int matrix_count = src->nhar * dst->vt_stride;
    for (int i = 0; i < matrix_count; ++i)
        dst->vt_magn[i] = src->vt_magn[i];

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
}

void lmix_frame_delete(lmix_frame *frame)
{
    if (frame == NULL)
        return;

    free(frame->vt_magn);
    free(frame->vs_har_phse);

    llsm_nosframe *noise = frame->noise;
    for (int i = 0; i < noise->nchannel; ++i) {
        llsm_sinframe *eenv = noise->eenv[i];

        if (eenv == NULL)
            continue;
        free(eenv->ampl);
        free(eenv->phse);
        free(eenv);
    }
    free(noise->eenv);
    free(noise->emin);
    free(noise->spec);
    free(noise);

    if (frame->wave != NULL) {
        free(frame->wave->samples);
        free(frame->wave);
    }
    free(frame);
}

void lmix_delete(lmix_data *lmix)
{
    if (lmix == NULL)
        return;

    free(lmix->conf.nosbandf);
    if (lmix->frames != NULL) {
        for (int i = 0; i < lmix->conf.nfrm; ++i)
            lmix_frame_delete(lmix->frames[i]);
        free(lmix->frames);
    }
    free(lmix);
}
