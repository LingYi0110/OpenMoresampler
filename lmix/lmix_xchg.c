#include "lmix.h"
#include "lmix_utils.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "../xchg.h"

// Convert LMIX parameters to an xchg node.
static xchg_node *lmix_params_to_xchg(llsm_conf conf)
{
    xchg_node *node = xchg_node_create(XCHG_TYPE_DICT, 11, 1);
    xchg_node *n_nfrm = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_nhop = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_nhar = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_nhare = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_nnos = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_noswarp = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_mvf = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_nosf = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_thop = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_nnosband = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_bands =
        xchg_node_create(XCHG_TYPE_FLOAT_ARRAY, conf.nnosband - 1, 1);

    n_nfrm->float_value = (float)conf.nfrm;
    n_nhop->float_value = (float)conf.nhop;
    n_nhar->float_value = (float)conf.nhar;
    n_nhare->float_value = (float)conf.nhare;
    n_nnos->float_value = (float)conf.nnos;
    n_noswarp->float_value = conf.noswarp;
    n_mvf->float_value = conf.mvf;
    n_nosf->float_value = conf.nosf;
    n_thop->float_value = conf.thop;
    n_nnosband->float_value = (float)conf.nnosband;

    int nband = n_bands->rows;
    if (nband > 0) {
        FP_TYPE *dst = (FP_TYPE *)n_bands->float_data;

        for (int i = 0; i != nband; ++i)
            dst[i] = conf.nosbandf[i];
    }

    xchg_node_set_field(node, n_nfrm, "_1", 0);
    xchg_node_set_field(node, n_nhop, "_2", 1);
    xchg_node_set_field(node, n_nhar, "_3", 2);
    xchg_node_set_field(node, n_nhare, "_4", 3);
    xchg_node_set_field(node, n_nnos, "_8", 4);
    xchg_node_set_field(node, n_noswarp, "_9", 5);
    xchg_node_set_field(node, n_mvf, "_a", 6);
    xchg_node_set_field(node, n_nosf, "_b", 7);
    xchg_node_set_field(node, n_thop, "_c", 8);
    xchg_node_set_field(node, n_nnosband, "_12", 9);
    xchg_node_set_field(node, n_bands, "_13", 10);
    return node;
}

static int16_t lmix_quant_phase(FP_TYPE phase)
{
    double p = phase;
    double turns = round(p * 0.1591549430918953);
    double wrapped = (p - turns * 6.283185307179586) * 10430.06004058427;

    return (int16_t)round(wrapped);
}

// Pack float deltas into a compact byte array.
static void lmix_pack_float_deltas(void *dst, const FP_TYPE *src, int n)
{
    unsigned char *out = (unsigned char *)dst;
    FP_TYPE first = src[0];
    if (n <= 1) {
        *(int *)out = 0;
        *(int *)(out + 4) = (int)0x80000000;
        *(FP_TYPE *)(out + 8) = first;
        return;
    }

    FP_TYPE min_d = 0.0f;
    FP_TYPE max_d = 0.0f;
    for (int i = 1; i < n; ++i) {
        FP_TYPE d = src[i] - src[i - 1];
        if (max_d < d)
            max_d = d;
        if (d < min_d)
            min_d = d;
    }
    FP_TYPE range = max_d - min_d;
    *(FP_TYPE *)out = min_d;
    *(FP_TYPE *)(out + 4) = range;
    *(FP_TYPE *)(out + 8) = first;
    FP_TYPE pred = first;
    FP_TYPE scale = 128.0f / range;
    FP_TYPE step = range * 0.0078125f;
    for (int i = 1; i < n; ++i) {
        int q = (int)round((src[i] - pred - min_d) * scale);
        if (q < -127)
            q = -127;
        if (q > 127)
            q = 127;
        out[12 + i - 1] = (unsigned char)q;
        pred = pred + (FP_TYPE)q * step + min_d;
    }
}

// Convert a noise frame to an xchg node.
static xchg_node *lmix_nosframe_to_xchg(llsm_nosframe *noise)
{
    int nch = noise->nchannel;
    int nnos = noise->nnos;
    int size = nnos + 11 + 2 * (3 * nch + 4);
    if (nch > 0) {
        for (int i = 0; i < nch; ++i)
            size += 4 * noise->eenv[i]->nhar;
    }

    xchg_node *node = xchg_node_create(XCHG_TYPE_DICT, 1, 1);
    xchg_node *blob = xchg_node_create(XCHG_TYPE_BYTE_ARRAY, size, 1);
    xchg_node_set_field(node, blob, "_40", 0);
    unsigned char *p = (unsigned char *)blob->byte_data;
    int *hdr = (int *)p;
    hdr[0] = nch;
    hdr[1] = nnos;
    int pack_off;
    if (nch <= 0) {
        pack_off = 8;
    } else {
        for (int i = 0; i < nch; ++i)
            hdr[i + 2] = noise->eenv[i]->nhar;
        pack_off = 4 * nch + 8;
    }
    lmix_pack_float_deltas(p + pack_off, noise->spec, nnos);
    int off = nnos + 11 + pack_off;
    if (nch > 0) {
        for (int i = 0; i < nch; ++i) {
            int16_t *w = (int16_t *)(p + off);
            w[0] = (int16_t)round(32767.0f * noise->emin[i]);
            off += 2;
        }
        for (int i = 0; i < nch; ++i) {
            llsm_sinframe *eenv = noise->eenv[i];
            int nh = eenv->nhar;

            if (nh <= 0)
                continue;
            int16_t *w = (int16_t *)(p + off);
            for (int j = 0; j < nh; ++j)
                w[j] = (int16_t)round(32767.0f * eenv->ampl[j]);
            off += 2 * nh;
            w = (int16_t *)(p + off);
            for (int j = 0; j < nh; ++j)
                w[j] = lmix_quant_phase(eenv->phse[j]);
            off += 2 * nh;
        }
    }
    return node;
}

// Convert a waveform to an xchg node.
static xchg_node *llsm_wave_to_xchg(llsm_wave *wave)
{
    xchg_node *node = xchg_node_create(XCHG_TYPE_DICT, 3, 1);
    xchg_node *n_fs = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_t0 = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    n_fs->float_value = (float)wave->fs;
    n_t0->float_value = wave->t0;
    xchg_node_set_field(node, n_fs, "_27", 0);
    xchg_node_set_field(node, n_t0, "_28", 1);
    xchg_node *blob =
        xchg_node_create(XCHG_TYPE_BYTE_ARRAY, 2 * wave->nsample + 14, 1);
    xchg_node_set_field(node, blob, "_40", 2);
    unsigned char *p = (unsigned char *)blob->byte_data;
    *(int *)p = wave->nsample;
    int n = wave->nsample;
    const FP_TYPE *src = wave->samples;
    FP_TYPE first = src[0];
    if (n <= 1) {
        *(int *)(p + 4) = 0;
        *(int *)(p + 8) = (int)0x80000000;
        *(FP_TYPE *)(p + 12) = first;
        return node;
    }
    FP_TYPE min_d = 0.0f;
    FP_TYPE max_d = 0.0f;
    for (int i = 1; i < n; ++i) {
        FP_TYPE d = src[i] - src[i - 1];
        if (max_d < d)
            max_d = d;
        if (d < min_d)
            min_d = d;
    }
    FP_TYPE range = max_d - min_d;
    *(FP_TYPE *)(p + 4) = min_d;
    *(FP_TYPE *)(p + 8) = range;
    *(FP_TYPE *)(p + 12) = first;
    FP_TYPE pred = first;
    FP_TYPE scale = 32768.0f / range;
    FP_TYPE step = range * 0.000030517578f;
    int16_t *w = (int16_t *)(p + 16);
    for (int i = 1; i < n; ++i) {
        int q = (int)round((src[i] - pred - min_d) * scale);
        if (q < -32767)
            q = -32767;
        if (q > 32767)
            q = 32767;
        w[i - 1] = (int16_t)q;
        pred = pred + (FP_TYPE)q * step + min_d;
    }
    return node;
}

// Convert a frame to an xchg node.
static xchg_node *lmix_frame_to_xchg(lmix_frame *frame)
{
    int has_wave = frame->wave != NULL;
    int nkeys = 2 + has_wave;
    if (frame->growl != 0.0f)
        nkeys += 1;
    if (frame->distortion != 0.0f)
        nkeys += 1;

    xchg_node *node = xchg_node_create(XCHG_TYPE_DICT, nkeys, 1);
    int nhar = frame->nhar;
    int blob_extra = 0;
    if (nhar > 0)
        blob_extra = nhar * frame->vt_stride + 11;
    xchg_node *blob = xchg_node_create(
        XCHG_TYPE_BYTE_ARRAY, blob_extra + 2 * nhar + 16, 1);
    unsigned char *p = (unsigned char *)blob->byte_data;
    *(FP_TYPE *)p = frame->f0;
    *(int16_t *)(p + 4) = (int16_t)nhar;
    *(int16_t *)(p + 6) = (int16_t)frame->vt_stride;
    *(FP_TYPE *)(p + 8) = frame->rd;
    *(FP_TYPE *)(p + 12) = frame->vt_freq_step;
    if (nhar > 0) {
        int16_t *ph = (int16_t *)(p + 16);
        for (int i = 0; i != nhar; ++i)
            ph[i] = lmix_quant_phase(frame->vs_har_phse[i]);
        lmix_pack_float_deltas(
            p + 16 + 2 * nhar, frame->vt_magn, nhar * frame->vt_stride);
    }
    xchg_node *nos = lmix_nosframe_to_xchg(frame->noise);
    xchg_node_set_field(node, blob, "_40", 0);
    xchg_node_set_field(node, nos, "_2f", 1);
    int idx_2b = 2;
    if (has_wave) {
        xchg_node_set_field(node, llsm_wave_to_xchg(frame->wave), "_5", 2);
        idx_2b = 3;
    }
    if (frame->growl != 0.0f) {
        xchg_node *growl_node = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
        growl_node->float_value = frame->growl;
        xchg_node_set_field(node, growl_node, "_2b", idx_2b);
        idx_2b += 1;
    }
    if (frame->distortion != 0.0f) {
        xchg_node *distortion_node = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
        distortion_node->float_value = frame->distortion;
        xchg_node_set_field(node, distortion_node, "_2d", idx_2b);
    }
    return node;
}

// Convert LMIX data to an xchg node.
xchg_node *lmix_to_xchg(lmix_data *lmix)
{
    xchg_node *root = xchg_node_create(XCHG_TYPE_DICT, 4, 1);
    xchg_node *n_version = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_duration = xchg_node_create(XCHG_TYPE_FLOAT, 0, 0);
    xchg_node *n_frames =
        xchg_node_create(XCHG_TYPE_OFFSET_ARRAY, lmix->conf.nfrm, 1);
    n_version->float_value = (float)lmix->version;
    n_duration->float_value = lmix->duration;

    int nfrm = n_frames->rows;
    if (nfrm > 0) {
        int i = 0;
        do {
            n_frames->children[i] =
                lmix_frame_to_xchg(lmix->frames[i]);
            ++i;
        } while (n_frames->rows > i);
    }

    xchg_node *n_params = lmix_params_to_xchg(lmix->conf);
    xchg_node_set_field(root, n_version, "version", 0);
    xchg_node_set_field(root, n_duration, "duration", 1);
    xchg_node_set_field(root, n_params, "_35", 2);
    xchg_node_set_field(root, n_frames, "_3a", 3);
    return root;
}

static void lmix_unpack_float_deltas(FP_TYPE *dst, const void *src, int n)
{
    const unsigned char *p = (const unsigned char *)src;
    FP_TYPE first = *(const FP_TYPE *)(p + 8);
    dst[0] = first;
    if (n <= 1)
        return;
    FP_TYPE min_d = *(const FP_TYPE *)p;
    FP_TYPE range = *(const FP_TYPE *)(p + 4);
    FP_TYPE step = range * 0.0078125f;
    FP_TYPE pred = first;
    for (int i = 1; i < n; ++i) {
        signed char q = (signed char)p[11 + i];
        pred = pred + (FP_TYPE)q * step + min_d;
        dst[i] = pred;
    }
}

// Read LMIX parameters from an xchg node.
static llsm_conf *lmix_params_from_xchg(xchg_node *node)
{
    if (xchg_node_has_fields(node, "_1", "_2", "_3", "_4", "_8", "_9",
                             "_a", "_b", "_c", "_12", "_13", NULL) == 0)
        return NULL;

    llsm_conf *conf = (llsm_conf *)calloc(1, sizeof(*conf));
    conf->nfrm = (int)xchg_node_get_field(node, "_1")->float_value;
    conf->nhop = (int)xchg_node_get_field(node, "_2")->float_value;
    conf->nhar = (int)xchg_node_get_field(node, "_3")->float_value;
    conf->nhare = (int)xchg_node_get_field(node, "_4")->float_value;
    conf->nnos = (int)xchg_node_get_field(node, "_8")->float_value;
    conf->noswarp = xchg_node_get_field(node, "_9")->float_value;
    conf->mvf = xchg_node_get_field(node, "_a")->float_value;
    conf->nosf = xchg_node_get_field(node, "_b")->float_value;
    conf->thop = xchg_node_get_field(node, "_c")->float_value;
    conf->nnosband = (int)xchg_node_get_field(node, "_12")->float_value;
    xchg_node *bands = xchg_node_get_field(node, "_13");
    int nband = bands->rows;
    conf->nosbandf = (FP_TYPE *)calloc((size_t)nband, sizeof(FP_TYPE));
    if (nband > 0) {
        FP_TYPE *src = (FP_TYPE *)bands->float_data;
        FP_TYPE *dst = conf->nosbandf;
        for (int i = 0; i < nband; ++i)
            dst[i] = src[i];
    }
    return conf;
}

// Read a sinusoidal frame from an xchg node.
static llsm_sinframe *lmix_sinframe_from_xchg(xchg_node *node)
{
    if (xchg_node_has_fields(node, "_40", NULL) != 0) {
        int *blob = (int *)xchg_node_get_field(node, "_40")->byte_data;

        int nhar = blob[0];
        llsm_sinframe *sinu = (llsm_sinframe *)malloc(sizeof(*sinu));
        sinu->nhar = nhar;
        sinu->ampl = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
        sinu->phse = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
        if (nhar > 0) {
            FP_TYPE *ampl_src = (FP_TYPE *)(blob + 1);
            FP_TYPE *phse_src = ampl_src + nhar;
            for (int i = 0; i != nhar; ++i)
                sinu->ampl[i] = ampl_src[i];
            for (int i = 0; i != nhar; ++i)
                sinu->phse[i] = phse_src[i];
        }
        return sinu;
    }
    if (xchg_node_has_fields(node, "_3", NULL) == 0)
        return NULL;
    int nhar = (int)xchg_node_get_field(node, "_3")->float_value;
    llsm_sinframe *sinu = (llsm_sinframe *)malloc(sizeof(*sinu));
    sinu->nhar = nhar;
    if (nhar <= 0) {
        sinu->ampl = NULL;
        sinu->phse = NULL;
        return sinu;
    }
    if (xchg_node_has_fields(node, "_14", "_15", NULL) == 0) {
        free(sinu);
        return NULL;
    }
    FP_TYPE *ampl_src = (FP_TYPE *)xchg_node_get_field(node, "_14")->float_data;
    FP_TYPE *phse_src = (FP_TYPE *)xchg_node_get_field(node, "_15")->float_data;
    sinu->ampl = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
    sinu->phse = (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
    for (int i = 0; i != nhar; ++i) {
        sinu->ampl[i] = ampl_src[i];
        sinu->phse[i] = phse_src[i];
    }
    return sinu;
}

// Read a noise frame from an xchg node.
static llsm_nosframe *lmix_nosframe_from_xchg(xchg_node *node, int version)
{
    if (xchg_node_has_fields(node, "_40", NULL) == 0) {
        if (xchg_node_has_fields(node, "_19", "_1d", "_1e", "_1f", "_8",
                                 NULL) == 0)
            return NULL;
        llsm_nosframe *noise = (llsm_nosframe *)malloc(sizeof(*noise));
        xchg_node *eenv_arr = xchg_node_get_field(node, "_19");
        xchg_node *emin_n = xchg_node_get_field(node, "_1d");
        xchg_node *spec_n = xchg_node_get_field(node, "_1e");
        int nch = (int)xchg_node_get_field(node, "_1f")->float_value;
        int nnos = (int)xchg_node_get_field(node, "_8")->float_value;
        noise->nchannel = nch;
        noise->nnos = nnos;
        int n_emin = emin_n->rows;
        int n_eenv = eenv_arr->rows;
        int n_spec = spec_n->rows;
        if (n_emin != nch || n_emin != n_eenv || nnos != n_spec) {
            free(noise);
            return NULL;
        }
        noise->spec = (FP_TYPE *)calloc((size_t)nnos, sizeof(FP_TYPE));
        noise->emin = (FP_TYPE *)calloc((size_t)nch, sizeof(FP_TYPE));
        FP_TYPE *spec_src = (FP_TYPE *)spec_n->float_data;
        FP_TYPE *emin_src = (FP_TYPE *)emin_n->float_data;
        for (int i = 0; i < nnos; ++i)
            noise->spec[i] = spec_src[i];
        for (int i = 0; i < nch; ++i)
            noise->emin[i] = emin_src[i];
        noise->eenv = (llsm_sinframe **)calloc(
            (size_t)nch, sizeof(llsm_sinframe *));
        for (int i = 0; i < nch; ++i)
            noise->eenv[i] = lmix_sinframe_from_xchg(
                (xchg_node *)eenv_arr->children[i]);
        return noise;
    }

    unsigned char *p =
        (unsigned char *)xchg_node_get_field(node, "_40")->byte_data;
    int *hdr = (int *)p;
    llsm_nosframe *noise = (llsm_nosframe *)malloc(sizeof(*noise));
    int nch = hdr[0];
    int nnos = hdr[1];
    noise->nchannel = nch;
    noise->nnos = nnos;
    noise->spec = (FP_TYPE *)malloc((size_t)nnos * sizeof(FP_TYPE));
    noise->emin = (FP_TYPE *)malloc((size_t)nch * sizeof(FP_TYPE));
    noise->eenv = (llsm_sinframe **)malloc((size_t)nch * sizeof(llsm_sinframe *));

    if (version > 500) {
        int pack_off;

        if (nch <= 0) {
            pack_off = 8;
        } else {
            for (int i = 0; i < nch; ++i) {
                int nh = hdr[i + 2];
                llsm_sinframe *eenv = (llsm_sinframe *)malloc(sizeof(*eenv));

                eenv->nhar = nh;
                noise->eenv[i] = eenv;
                if (nh > 0) {
                    eenv->ampl = (FP_TYPE *)malloc((size_t)nh * sizeof(FP_TYPE));
                    eenv->phse = (FP_TYPE *)malloc((size_t)nh * sizeof(FP_TYPE));
                } else {
                    eenv->ampl = NULL;
                    eenv->phse = NULL;
                }
            }
            pack_off = 4 * nch + 8;
        }
        lmix_unpack_float_deltas(noise->spec, p + pack_off, nnos);
        int off = pack_off + nnos + 11;
        if (nch > 0) {
            for (int i = 0; i < nch; ++i) {
                int16_t *w = (int16_t *)(p + off);
                noise->emin[i] = (FP_TYPE)w[0] * 0.000030518509f;
                off += 2;
            }
            for (int i = 0; i < nch; ++i) {
                llsm_sinframe *eenv = noise->eenv[i];
                int nh = eenv->nhar;

                if (nh <= 0)
                    continue;
                int16_t *w = (int16_t *)(p + off);
                for (int j = 0; j < nh; ++j)
                    eenv->ampl[j] = (FP_TYPE)w[j] * 0.000030518509f;
                off += 2 * nh;
                w = (int16_t *)(p + off);
                for (int j = 0; j < nh; ++j)
                    eenv->phse[j] =
                        (FP_TYPE)w[j] * 0.000030518509f * 3.141592653589793f;
                off += 2 * nh;
            }
        }
        return noise;
    }

    int off = nch > 0 ? 4 * nch + 12 : 8;
    if (nch > 0) {
        for (int i = 0; i < nch; ++i) {
            int nh = hdr[i + 2];
            llsm_sinframe *eenv = (llsm_sinframe *)malloc(sizeof(*eenv));

            eenv->nhar = nh;
            noise->eenv[i] = eenv;
            if (nh > 0) {
                eenv->ampl = (FP_TYPE *)malloc((size_t)nh * sizeof(FP_TYPE));
                eenv->phse = (FP_TYPE *)malloc((size_t)nh * sizeof(FP_TYPE));
            } else {
                eenv->ampl = NULL;
                eenv->phse = NULL;
            }
        }
    }
    if (nnos > 0) {
        FP_TYPE *src = (FP_TYPE *)(p + off);
        for (int i = 0; i < nnos; ++i)
            noise->spec[i] = src[i];
        off += 4 * nnos;
    }
    if (nch > 0) {
        FP_TYPE *src = (FP_TYPE *)(p + off);
        for (int i = 0; i < nch; ++i)
            noise->emin[i] = src[i];
        off += 4 * nch;
        for (int i = 0; i < nch; ++i) {
            llsm_sinframe *eenv = noise->eenv[i];
            int nh = eenv->nhar;
            FP_TYPE *srcf;

            if (nh <= 0)
                continue;
            srcf = (FP_TYPE *)(p + off);
            for (int j = 0; j < nh; ++j)
                eenv->ampl[j] = srcf[j];
            off += 4 * nh;
            srcf = (FP_TYPE *)(p + off);
            for (int j = 0; j < nh; ++j)
                eenv->phse[j] = srcf[j];
            off += 4 * nh;
        }
    }
    return noise;
}

// Read a waveform from an xchg node.
static llsm_wave *llsm_wave_from_xchg(xchg_node *node)
{
    llsm_wave *wave = (llsm_wave *)malloc(sizeof(*wave));
    wave->fs = (int)xchg_node_get_field(node, "_27")->float_value;
    wave->t0 = xchg_node_get_field(node, "_28")->float_value;
    if (xchg_node_has_fields(node, "_40", NULL) == 0) {
        xchg_node *arr = xchg_node_get_field(node, "_11");

        int n = arr->rows;
        wave->nsample = n;
        wave->samples = (FP_TYPE *)calloc((size_t)n, sizeof(FP_TYPE));
        FP_TYPE *src = (FP_TYPE *)arr->float_data;
        for (int i = 0; i < n; ++i)
            wave->samples[i] = src[i];
        return wave;
    }
    {
        unsigned char *p =
            (unsigned char *)xchg_node_get_field(node, "_40")->byte_data;

        int n = *(int *)p;
        wave->nsample = n;
        FP_TYPE *dst = (FP_TYPE *)calloc((size_t)n, sizeof(FP_TYPE));
        wave->samples = dst;
        FP_TYPE min_d = *(FP_TYPE *)(p + 4);
        FP_TYPE range = *(FP_TYPE *)(p + 8);
        FP_TYPE pred = *(FP_TYPE *)(p + 12);
        dst[0] = pred;
        if (n <= 1)
            return wave;
        FP_TYPE step = range * 0.000030517578f;
        int16_t *w = (int16_t *)(p + 16);
        for (int i = 1; i < n; ++i) {
            pred = pred + (FP_TYPE)w[i - 1] * step + min_d;
            dst[i] = pred;
        }
        return wave;
    }
}

static void lmix_frame_free_partial(lmix_frame *frame)
{
    if (frame->nhar > 0) {
        free(frame->vs_har_phse);
        free(frame->vt_magn);
    }
    free(frame);
}

// Read a frame from an xchg node.
static lmix_frame *lmix_frame_from_xchg(xchg_node *node, int version)
{
    if (xchg_node_has_fields(node, "_2f", NULL) == 0)
        return NULL;
    lmix_frame *frame = (lmix_frame *)calloc(1, sizeof(*frame));
    if (xchg_node_has_fields(node, "_40", NULL) != 0) {
        unsigned char *p =
            (unsigned char *)xchg_node_get_field(node, "_40")->byte_data;
        frame->f0 = *(FP_TYPE *)p;
        if (version > 500) {
            int nhar = *(int16_t *)(p + 4);
            int stride = *(int16_t *)(p + 6);
            frame->rd = *(FP_TYPE *)(p + 8);
            frame->vt_freq_step = *(FP_TYPE *)(p + 12);
            frame->nhar = nhar;
            frame->vt_stride = stride;
            if (nhar > 0) {
                int16_t *ph = (int16_t *)(p + 16);

                frame->vs_har_phse =
                    (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
                int nmagn = nhar * stride;
                frame->vt_magn =
                    (FP_TYPE *)calloc((size_t)nmagn, sizeof(FP_TYPE));
                for (int i = 0; i != nhar; ++i)
                    frame->vs_har_phse[i] =
                        (FP_TYPE)ph[i] * 0.000030518509f * 3.141592653589793f;
                lmix_unpack_float_deltas(
                    frame->vt_magn, p + 16 + 2 * nhar, nmagn);
            }
        } else {
            int nhar = *(int *)(p + 4);
            int stride = *(int *)(p + 8);
            frame->rd = *(FP_TYPE *)(p + 12);
            frame->vt_freq_step = *(FP_TYPE *)(p + 16);
            frame->nhar = nhar;
            frame->vt_stride = stride;
            if (nhar > 0) {
                FP_TYPE *ph = (FP_TYPE *)(p + 20);

                frame->vs_har_phse =
                    (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
                int nmagn = nhar * stride;
                frame->vt_magn =
                    (FP_TYPE *)calloc((size_t)nmagn, sizeof(FP_TYPE));
                for (int i = 0; i != nhar; ++i)
                    frame->vs_har_phse[i] = ph[i];
                for (int i = 0; i < nmagn; ++i)
                    frame->vt_magn[i] =
                        ((FP_TYPE *)(p + 4 * nhar + 20))[i];
            }
        }
    } else {
        if (xchg_node_has_fields(node, "_20", "_3", "_25", "_2a", "_2c",
                                 NULL) == 0) {
            free(frame);
            return NULL;
        }
        frame->f0 = xchg_node_get_field(node, "_20")->float_value;
        int nhar = (int)xchg_node_get_field(node, "_3")->float_value;
        int stride = (int)xchg_node_get_field(node, "_25")->float_value;
        frame->rd = xchg_node_get_field(node, "_2a")->float_value;
        frame->vt_freq_step = xchg_node_get_field(node, "_2c")->float_value;
        frame->nhar = nhar;
        frame->vt_stride = stride;
        if (nhar > 0) {
            if (xchg_node_has_fields(node, "_31", "_32", NULL) == 0) {
                free(frame);
                return NULL;
            }
            xchg_node *ph_n = xchg_node_get_field(node, "_31");
            xchg_node *mg_n = xchg_node_get_field(node, "_32");
            int nph = ph_n->rows;
            int nmagn = mg_n->rows;
            if (nph != nhar || nmagn != nhar * stride) {
                free(frame);
                return NULL;
            }
            frame->vs_har_phse =
                (FP_TYPE *)calloc((size_t)nhar, sizeof(FP_TYPE));
            frame->vt_magn =
                (FP_TYPE *)calloc((size_t)nmagn, sizeof(FP_TYPE));
            FP_TYPE *src = (FP_TYPE *)ph_n->float_data;
            for (int i = 0; i < nhar; ++i)
                frame->vs_har_phse[i] = src[i];
            src = (FP_TYPE *)mg_n->float_data;
            for (int i = 0; i < nmagn; ++i)
                frame->vt_magn[i] = src[i];
        }
    }

    xchg_node *wave_n = xchg_node_get_field(node, "_5");
    if (wave_n != NULL) {
        frame->wave = llsm_wave_from_xchg(wave_n);
        if (frame->wave == NULL) {
            lmix_frame_free_partial(frame);
            return NULL;
        }
    }
    frame->noise = lmix_nosframe_from_xchg(
        xchg_node_get_field(node, "_2f"), version);
    if (frame->noise == NULL) {
        if (frame->wave != NULL) {
            free(frame->wave->samples);
            free(frame->wave);
        }
        lmix_frame_free_partial(frame);
        return NULL;
    }
    xchg_node *growl_node = xchg_node_get_field(node, "_2b");
    xchg_node *distortion_node = xchg_node_get_field(node, "_2d");
    if (growl_node != NULL)
        frame->growl = growl_node->float_value;
    if (distortion_node != NULL)
        frame->distortion = distortion_node->float_value;
    return frame;
}

// Read LMIX data from an xchg node.
lmix_data *lmix_from_xchg(xchg_node *node)
{
    if (xchg_node_has_fields(node, "version", "duration", "_35", "_3a",
                             NULL) == 0)
        return NULL;
    lmix_data *lmix = (lmix_data *)malloc(sizeof(*lmix));
    int version = (int)xchg_node_get_field(node, "version")->float_value;
    lmix->version = version;
    lmix->duration = xchg_node_get_field(node, "duration")->float_value;
    llsm_conf *conf = lmix_params_from_xchg(
        xchg_node_get_field(node, "_35"));
    if (conf == NULL) {
        free(lmix);
        return NULL;
    }
    lmix->conf = *conf;
    free(conf);
    xchg_node *frames_n = xchg_node_get_field(node, "_3a");
    int nfrm = frames_n->rows;
    if (nfrm != lmix->conf.nfrm) {
        free(lmix->conf.nosbandf);
        free(lmix);
        return NULL;
    }
    lmix->frames = (lmix_frame **)calloc((size_t)nfrm, sizeof(lmix_frame *));
    for (int i = 0; i < nfrm; ++i) {
        xchg_node *child = (xchg_node *)frames_n->children[i];

        if (child == NULL) {
            lmix->frames[i] = NULL;
            continue;
        }
        lmix->frames[i] = lmix_frame_from_xchg(child, version);
        if (lmix->frames[i] == NULL) {
            for (int j = 0; j < i; ++j)
                lmix_frame_delete(lmix->frames[j]);
            free(lmix->frames);
            free(lmix->conf.nosbandf);
            free(lmix);
            return NULL;
        }
    }
    return lmix;
}
