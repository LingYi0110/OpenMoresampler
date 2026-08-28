#include "otogen.h"

#include <math.h>

#include "../logging.h"
#include "../utils.h"
#include "../lmix/lmix_utils.h"

otogen_ctx *otogen_ctx_create(void)
{
    otogen_ctx *ctx = (otogen_ctx *)calloc(1, sizeof(otogen_ctx));
    lrh_inference_stprune = 10.0f;
    ctx->nstate = 3;
    ctx->nstream = 1;
    ctx->ncep = 12;
    ctx->nmel = 20;
    ctx->nmel_delta = 18;
    ctx->fs = 44000.0f;
    ctx->fmin = 50.0f;
    ctx->fmax = 8000.0f;
    ctx->hop_sec = 0.01f;
    ctx->win_sec = 0.025f;
    ctx->fft_sec = 0.4f;
    ctx->silence_thresh = 0.005f;
    ctx->pad_sec = 0.1f;
    ctx->silence_floor = 0.5f;
    ctx->nskip = 10;
    return ctx;
}

void otogen_ctx_free(otogen_ctx *ctx)
{
    if (ctx == NULL)
        return;
    otogen_g2p_free(ctx->g2p);
    otogen_phdict_free(ctx->phdict);
    if (ctx->wav_paths != NULL) {
        for (int i = 0; i < ctx->nfile; i++)
            free(ctx->wav_paths[i]);
        free(ctx->wav_paths);
    }
    if (ctx->aliases != NULL) {
        for (int i = 0; i < ctx->nfile; i++)
            free(ctx->aliases[i]);
        free(ctx->aliases);
    }
    if (ctx->labels != NULL) {
        for (int i = 0; i < ctx->nfile; i++)
            otogen_label_free(ctx->labels[i]);
        free(ctx->labels);
    }
    free(ctx->infer_nt);
    free(ctx->trim);
    if (ctx->band_e != NULL) {
        for (int i = 0; i < ctx->nfile; i++) {
            int nt = (ctx->observset != NULL && ctx->observset->samples[i] != NULL)
                ? ctx->observset->samples[i]->nt : 0;
            if (ctx->band_e[i] == NULL)
                continue;
            for (int t = 0; t < nt; t++)
                free(ctx->band_e[i][t]);
            free(ctx->band_e[i]);
        }
        free(ctx->band_e);
    }
    lrh_delete_observset(ctx->observset);
    lrh_delete_segset(ctx->segset);
    lrh_delete_model(ctx->model);
    otogen_bmodel_free(ctx->boundary_model);
    free(ctx);
}

int otogen_scan_wavs(otogen_ctx *ctx, const wchar_t *dir)
{
    size_t n = 0;
    wchar_t **names = list_files(dir, L".wav", &n);
    if (names == NULL)
        return 0;
    wchar_t *pref = path_dir_prefix(dir);
    ctx->wav_paths = (wchar_t **)realloc(ctx->wav_paths,
        (size_t)(ctx->nfile + (int)n) * sizeof(wchar_t *));
    for (int i = 0; i < (int)n; i++) {
        ctx->wav_paths[ctx->nfile + i] = wstr_concat(pref, names[i]);
        free(names[i]);
    }
    ctx->nfile += (int)n;
    free(names);
    free(pref);
    return 1;
}

int otogen_aliases_from_wavs(otogen_ctx *ctx)
{
    if (ctx->wav_paths == NULL) {
        fwrite("Error: list of file names not initialized yet.\n", 1, 0x2F, stderr);
        return 0;
    }
    if (ctx->aliases != NULL) {
        for (int i = 0; i < ctx->nfile; i++)
            free(ctx->aliases[i]);
        free(ctx->aliases);
    }
    ctx->aliases = (wchar_t **)calloc((size_t)ctx->nfile, sizeof(wchar_t *));
    for (int i = 0; i < ctx->nfile; i++) {
        if (ctx->wav_paths[i] == NULL)
            return 0;
        wchar_t *base = path_basename(ctx->wav_paths[i]);
        size_t n = wcslen(base);
        ctx->aliases[i] = wstr_ndup(base, n >= 4 ? n - 4 : n);
        free(base);
    }
    return 1;
}

static int pow2_round(FP_TYPE x)
{
    return (int)pow(2.0, round(log2((double)x)));
}

static FP_TYPE **stft_hamming(FP_TYPE *x, int nx, int fs, FP_TYPE hop_sec, FP_TYPE fft_sec, int *nfrm, int *nfft)
{
    FP_TYPE hop_f = (FP_TYPE)fs * hop_sec;
    int fft = pow2_round(fft_sec * (FP_TYPE)fs);
    int nwin_i = (int)ceil((double)((FP_TYPE)fs * fft_sec));

    *nfrm = (int)round((double)((FP_TYPE)nx / hop_f));
    *nfft = fft;
    FP_TYPE **magn = (FP_TYPE **)calloc((size_t)*nfrm, sizeof(FP_TYPE *));
    for (int i = 0; i < *nfrm; i++)
        magn[i] = (FP_TYPE *)calloc((size_t)((fft >> 1) + 1), sizeof(FP_TYPE));
    int *center = (int *)calloc((size_t)*nfrm, sizeof(int));
    int *nwin = (int *)calloc((size_t)*nfrm, sizeof(int));
    for (int i = 0; i < *nfrm; i++) {
        center[i] = (int)round((double)((FP_TYPE)i * hop_f));
        nwin[i] = nwin_i;
    }
    cig_stft_forward(x, nx, center, nwin, *nfrm, fft, "hamming", 1, 2, NULL, NULL, magn, NULL);
    free(center);
    free(nwin);
    return magn;
}

static void free_spgm(FP_TYPE **S, int nfrm)
{
    if (S == NULL)
        return;
    for (int i = 0; i < nfrm; i++)
        free(S[i]);
    free(S);
}

static void trim_silence(FP_TYPE *x, int *nx, FP_TYPE fs, FP_TYPE thresh, FP_TYPE *t0, FP_TYPE *t1)
{
    int n = *nx;
    int i0 = 0;
    int trail = 0;
    FP_TYPE start = 0.0f;
    FP_TYPE endt = 0.0f;
    int i;

    if (n > 0 && x[0] < thresh) {
        i = 0;
        while (++i != n) {
            if (x[i] >= thresh) {
                i0 = i;
                start = (FP_TYPE)i;
                goto found_head;
            }
        }
        i0 = 0;
        start = 0.0f;
    }
found_head:
    i = n - 1;
    if (n - 1 >= 0) {
        do {
            if (x[i] >= thresh) {
                trail = n - i;
                endt = (FP_TYPE)trail;
                goto found_tail;
            }
            --i;
        } while (i != -1);
    }
    trail = 0;
    endt = 0.0f;
found_tail:
    *t0 = start / fs;
    *t1 = endt / fs;
    if (n - i0 > 0) {
        FP_TYPE *dst = x;
        do {
            *dst = dst[i0];
            ++dst;
        } while (dst != x + (n - i0));
    }
    *nx = n - (trail + i0);
}

int otogen_analyze_audio(otogen_ctx *ctx)
{
    if (ctx->wav_paths == NULL) {
        fwrite("Error: list of file names not initialized yet.\n", 1, 0x2F, stderr);
        return 0;
    }
    if (ctx->observset != NULL)
        lrh_delete_observset(ctx->observset);

    int nfft_fb = pow2_round(ctx->fft_sec * ctx->fs);
    filterbank *fbank = cig_create_melfreq_filterbank(nfft_fb / 2 + 1, ctx->fs * 0.5f,
        ctx->nmel, ctx->fmin, ctx->fmax, 1.0f, 400.0f);
    filterbank *fbank_d = cig_create_melfreq_filterbank(nfft_fb / 2 + 1, ctx->fs * 0.5f,
        ctx->nmel_delta, ctx->fmin, ctx->fmax, 1.0f, 400.0f);

    ctx->observset = lrh_create_empty_observset(ctx->nfile);
    ctx->band_e = (FP_TYPE ***)calloc((size_t)ctx->nfile, sizeof(FP_TYPE **));
    ctx->trim = (FP_TYPE *)calloc((size_t)ctx->nfile, 2 * sizeof(FP_TYPE));

    int nstream = ctx->nstream + 1;
    int ndim = ctx->ncep + 1;

    for (int i = 0; i < ctx->nfile; i++) {
        char *path_mbs;
        FP_TYPE *wav;
        int fs = 0, nbit = 0, nx = 0;
        FP_TYPE *rs;
        int nrs = 0;
        int pad;
        FP_TYPE *padded;
        int npad;
        FP_TYPE peak;
        FP_TYPE **spec;
        int nfrm = 0, nfft = 0;
        FP_TYPE **be;
        FP_TYPE **be18;
        FP_TYPE **cc;
        int *ndims;
        lrh_observ *ob;

        path_mbs = wchar_to_mbs(ctx->wav_paths[i]);
        wav = wavread(path_mbs, &fs, &nbit, &nx);
        if (wav == NULL) {
            fprintf(stderr, "Error: can't open wav file '%s'.\n", path_mbs);
            free(path_mbs);
            lrh_delete_observset(ctx->observset);
            cig_delete_filterbank(fbank);
            return 0;
        }
        free(path_mbs);
        trim_silence(wav, &nx, (FP_TYPE)fs, ctx->silence_thresh,
            &ctx->trim[2 * i], &ctx->trim[2 * i + 1]);
        printf("\rAnalyzing %d/%d...", i + 1, ctx->nfile);
        fflush(stdout);

        rs = cig_rresample(wav, nx, ctx->fs / (FP_TYPE)fs, &nrs);
        pad = (int)round((double)(ctx->pad_sec * ctx->fs));
        npad = 2 * pad + nrs;
        padded = (FP_TYPE *)calloc((size_t)npad, sizeof(FP_TYPE));
        peak = 0.0f;
        if (nrs > 0) {
            for (int t = 0; t < nrs; t++) {
                FP_TYPE a = (FP_TYPE)fabs((double)rs[t]);
                if (peak < a)
                    peak = a;
            }
            for (int t = 0; t < nrs; t++)
                padded[pad + t] = rs[t] * (1.0f / peak);
        }
        for (int t = 0; t < npad; t++) {
            int r = rand();
            padded[t] += ((FP_TYPE)r * 0.000030518509f - 0.5f) * 2.0f * 0.005f;
        }

        spec = stft_hamming(padded, npad, (int)ctx->fs, ctx->hop_sec, ctx->fft_sec, &nfrm, &nfft);
        be = cig_filterbank_spectrogram(fbank, spec, nfrm, nfft, (int)ctx->fs, 1);
        be18 = cig_filterbank_spectrogram(fbank_d, spec, nfrm, nfft, (int)ctx->fs, 1);
        cc = (FP_TYPE **)malloc((size_t)nfrm * sizeof(FP_TYPE *));
        {
            int ncc = ctx->ncep < ctx->nmel - 1 ? ctx->ncep : ctx->nmel - 1;
            for (int j = 0; j < nfrm; j++) {
                FP_TYPE *dct = cig_dct(be[j], ctx->nmel);
                FP_TYPE dc = dct[0];
                if (ncc > 0)
                    memmove(dct, dct + 1, (size_t)ncc * sizeof(FP_TYPE));
                dct[ncc] = dc;
                cc[j] = dct;
            }
        }
        ctx->band_e[i] = be18;

        ndims = (int *)malloc((size_t)nstream * sizeof(int));
        for (int s = 0; s < nstream; s++)
            ndims[s] = ndim;
        ob = lrh_create_observ(nstream, nfrm, ndims);
        free(ndims);

        for (int j = 0; j < nfrm; j++) {
            for (int k = 0; k < ctx->ncep; k++)
                lrh_obm(ob, j, k, 0) = cc[j][k];
            lrh_obm(ob, j, ctx->ncep, 0) = cc[j][ctx->ncep];
        }
        for (int s = 1; s < nstream; s++) {
            for (int j = 0; j < nfrm; j++) {
                for (int k = 0; k < ndim; k++) {
                    if (j == 0 || j == nfrm - 1)
                        lrh_obm(ob, j, k, s) = 0.0f;
                    else
                        lrh_obm(ob, j, k, s) =
                            0.5f * (lrh_obm(ob, j + 1, k, s - 1) - lrh_obm(ob, j - 1, k, s - 1));
                }
            }
        }

        ctx->observset->samples[i] = ob;

        free_spgm(spec, nfrm);
        free_spgm(be, nfrm);
        free_spgm(cc, nfrm);
        free(wav);
        free(rs);
        free(padded);
    }
    putchar('\n');
    cig_delete_filterbank(fbank);
    cig_delete_filterbank(fbank_d);

    {
        int ns = ctx->nstream + 1;
        int nd = ctx->ncep + 1;

        for (int d = 0; d < ns; d++) {
            FP_TYPE *sum = (FP_TYPE *)calloc((size_t)nd, sizeof(FP_TYPE));
            FP_TYPE *sq = (FP_TYPE *)calloc((size_t)nd, sizeof(FP_TYPE));
            int ntot = 0;

            for (int f = 0; f < ctx->nfile; f++) {
                lrh_observ *ob = ctx->observset->samples[f];
                int ntd = ob->ndim[d];
                ntot += ob->nt;
                if (ob->nt <= 0 || ntd <= 0)
                    continue;
                for (int t = 0; t < ob->nt; t++) {
                    for (int k = 0; k < ntd; k++) {
                        FP_TYPE v = lrh_obm(ob, t, k, d);
                        sum[k] += v;
                        sq[k] += v * v;
                    }
                }
            }
            {
                FP_TYPE inv = 1.0f / (FP_TYPE)ntot;
                for (int k = 0; k < nd; k++) {
                    FP_TYPE mean = sum[k] * inv;
                    sum[k] = mean;
                    sq[k] = (FP_TYPE)sqrt((double)(sq[k] * inv - mean * mean));
                }
            }
            for (int f = 0; f < ctx->nfile; f++) {
                lrh_observ *ob = ctx->observset->samples[f];
                if (ob->nt <= 0)
                    continue;
                for (int t = 0; t < ob->nt; t++) {
                    for (int k = 0; k < nd; k++) {
                        FP_TYPE v = lrh_obm(ob, t, k, d) - sum[k];
                        lrh_obm(ob, t, k, d) = v / sq[k] * 10.0f;
                    }
                }
            }
            free(sum);
            free(sq);
        }
    }
    return 1;
}

void otogen_init_segmentation(otogen_ctx *ctx)
{
    int nstream;

    if (ctx->observset == NULL) {
        fwrite("Error: audio data not initialized yet.\n", 1, 0x27, stderr);
        return;
    }
    if (ctx->phdict == NULL) {
        fwrite("Error: phoneme dictionary not initialized yet.\n", 1, 0x2F, stderr);
        return;
    }
    if (ctx->labels == NULL) {
        if (ctx->g2p == NULL) {
            fwrite("Error: grapheme-to-phoneme dictionary not initialized yet.\n", 1, 0x3B, stderr);
            return;
        }
        if (ctx->aliases == NULL) {
            fwrite("Error: list of plain text not initialized yet.\n", 1, 0x2F, stderr);
            return;
        }
        ctx->labels = (otogen_label **)calloc((size_t)ctx->nfile, sizeof(otogen_label *));
        for (int i = 0; i < ctx->nfile; i++)
            ctx->labels[i] = otogen_g2p_apply(ctx->g2p, ctx->aliases[i]);
    }

    if (ctx->segset != NULL)
        lrh_delete_segset(ctx->segset);
    ctx->segset = lrh_create_empty_segset(ctx->nfile);
    ctx->infer_nt = (int *)calloc((size_t)ctx->nfile, sizeof(int));
    nstream = ctx->nstream + 1;

    for (int i = 0; i < ctx->nfile; i++) {
        otogen_label *lab = ctx->labels[i];
        lrh_observ *ob = ctx->observset->samples[i];
        int nph = lab->nseg;
        int nseg = nph * ctx->nstate;
        FP_TYPE duration = (FP_TYPE)ob->nt * ctx->hop_sec;
        int nnonconsonant = 0;
        lrh_seg *seg;

        for (int p = 0; p < nph; p++) {
            lab->tbegin[p] = (FP_TYPE)p * duration / (FP_TYPE)nph;
            lab->tend[p] = (FP_TYPE)(p + 1) * duration / (FP_TYPE)nph;
            nnonconsonant += !otogen_is_consonant(ctx->phdict, lab->phoneme[p]);
        }

        if (nnonconsonant - 2 <= 0) {
            ctx->infer_nt[i] = 50;
        } else {
            ctx->infer_nt[i] = (int)((duration - (ctx->fft_sec + ctx->fft_sec))
                / ((FP_TYPE)ctx->nstate * ctx->hop_sec
                    * (FP_TYPE)(nnonconsonant - 2)));
        }
        seg = lrh_create_seg(nstream, nseg);
        for (int p = 0; p < nph; p++) {
            otogen_hnode *hn = otogen_hash_find(ctx->phdict->all, lab->phoneme[p]);
            int base;

            if (hn == NULL) {
                char *m = wchar_to_mbs(lab->phoneme[p]);
                fprintf(stderr, "Error: unknown phoneme '%s' at %d.\n", m, i);
                free(m);
                return;
            }
            base = hn->index * ctx->nstate;
            for (int st = 0; st < ctx->nstate; st++) {
                int si = p * ctx->nstate + st;
                seg->time[si] = (int)floor(
                    (lab->tbegin[p]
                        + (lab->tend[p] - lab->tbegin[p])
                            / (FP_TYPE)ctx->nstate * (FP_TYPE)(st + 1))
                    / ctx->hop_sec);
                seg->durstate[si] = base + st;
                for (int s = 0; s < nstream; s++)
                    seg->outstate[s][si] = base + st;
            }
        }
        ctx->segset->samples[i] = seg;
    }
}

void otogen_init_model(otogen_ctx *ctx, int skip_stats)
{
    int nstream;
    int ndur;
    int *ngmm;
    int *nmix;
    int *ndim;
    int feat_dim;
    lrh_model_stat *st;

    if (ctx->observset == NULL) {
        fwrite("Error: audio data not initialized yet.\n", 1, 0x27, stderr);
        return;
    }
    if (ctx->phdict == NULL) {
        fwrite("Error: phoneme dictionary not initialized yet.\n", 1, 0x2F, stderr);
        return;
    }
    if (ctx->segset == NULL) {
        fwrite("Error: segmentation not initialized yet.\n", 1, 0x29, stderr);
        return;
    }
    if (ctx->model != NULL)
        lrh_delete_model(ctx->model);

    nstream = ctx->nstream + 1;
    ndur = ctx->nstate * ctx->phdict->all->count;
    feat_dim = ctx->ncep + 1;
    ngmm = (int *)malloc((size_t)nstream * sizeof(int));
    nmix = (int *)malloc((size_t)nstream * sizeof(int));
    ndim = (int *)malloc((size_t)nstream * sizeof(int));
    for (int i = 0; i < nstream; i++) {
        ngmm[i] = ndur;
        nmix[i] = 1;
        ndim[i] = feat_dim;
    }
    ctx->model = lrh_create_model(nstream, ngmm, ndur, nmix, ndim);
    free(ngmm);
    free(nmix);
    free(ndim);

    for (int i = 0; i < ctx->phdict->all->count; i++) {
        if (otogen_is_silence(ctx->phdict, ctx->phdict->names[i])) {
            for (int stt = 0; stt < ctx->nstate; stt++) {
                int di = i * ctx->nstate + stt;
                ctx->model->durations[di]->vfloor = 400.0f;
            }
        }
    }

    st = lrh_model_stat_from_model(ctx->model);
    for (int i = 0; i < ctx->nfile; i++)
        lrh_collect_init(st, ctx->observset->samples[i], ctx->segset->samples[i]);
    lrh_model_update(ctx->model, st, 0);

    if (skip_stats == 0) {
        int nfile = ctx->nfile;
        struct otogen_stream_acc {
            FP_TYPE *valsum;
            FP_TYPE *sqrsum;
            int ntot;
            int ndim;
        } **acc;

        acc = (struct otogen_stream_acc **)calloc((size_t)nstream, sizeof(*acc));
        for (int s = 0; s < nstream; s++) {
            acc[s] = (struct otogen_stream_acc *)malloc(sizeof(struct otogen_stream_acc));
            acc[s]->valsum = (FP_TYPE *)calloc((size_t)feat_dim, sizeof(FP_TYPE));
            acc[s]->sqrsum = (FP_TYPE *)calloc((size_t)feat_dim, sizeof(FP_TYPE));
            acc[s]->ntot = 0;
            acc[s]->ndim = feat_dim;
            for (int f = 0; f < nfile; f++) {
                lrh_observ *ob = ctx->observset->samples[f];
                int nd;

                if (ob->nt > 0) {
                    nd = ob->ndim[s];
                    if (nd > 0) {
                        for (int t = 0; t < ob->nt; t++) {
                            for (int k = 0; k < nd; k++) {
                                FP_TYPE v = lrh_obm(ob, t, k, s);
                                acc[s]->valsum[k] += v;
                                acc[s]->sqrsum[k] += v * v;
                            }
                        }
                    }
                }
                acc[s]->ntot += ob->nt;
            }
        }
        if (ndur > 0) {
            for (int s = 0; s < nstream; s++) {
                FP_TYPE inv;
                lrh_stream *stm;

                if (feat_dim <= 0)
                    continue;
                inv = 1.0f / (FP_TYPE)acc[s]->ntot;
                stm = ctx->model->streams[s];
                for (int g = 0; g < ndur; g++) {
                    lrh_gmm *gmm = stm->gmms[g];
                    for (int k = 0; k < feat_dim; k++) {
                        FP_TYPE mean = acc[s]->valsum[k] * inv;
                        FP_TYPE m2 = acc[s]->sqrsum[k] * inv;
                        gmm->mean[k] = mean;
                        gmm->var[k] = m2 - mean * mean;
                    }
                }
            }
        }
        for (int s = 0; s < nstream; s++) {
            if (acc[s] != NULL) {
                free(acc[s]->valsum);
                free(acc[s]->sqrsum);
                free(acc[s]);
            }
        }
        free(acc);
    }
    lrh_delete_model_stat(st);
}

void otogen_pool_vowel_gmms(otogen_ctx *ctx)
{
    lrh_model *model = ctx->model;
    int nstream = model->nstream;
    int nstate = ctx->nstate;
    int ngmm;
    lrh_gmm_stat **gmm_stats;
    lrh_duration_stat **dur_stats;

    if (nstream <= 0)
        return;
    ngmm = model->streams[0]->ngmm;
    gmm_stats = (lrh_gmm_stat **)malloc((size_t)nstream * sizeof(lrh_gmm_stat *));
    dur_stats = (lrh_duration_stat **)malloc((size_t)nstate * sizeof(lrh_duration_stat *));
    for (int s = 0; s < nstream; s++)
        gmm_stats[s] = lrh_gmm_stat_from_gmm(model->streams[s]->gmms[0]);
    for (int i = 0; i < nstate; i++)
        dur_stats[i] = lrh_duration_stat_from_duration(model->durations[0]);

    if (ngmm > 0) {
        for (int s = 0; s < nstream; s++) {
            lrh_stream *stm = model->streams[s];
            for (int i = 0; i < ngmm; i++) {
                int ph = i / nstate;
                int stt = i % nstate;
                wchar_t *name = ctx->phdict->names[ph];
                lrh_duration_stat *ds;
                lrh_gmm_stat *gs;
                lrh_gmm *gmm;
                int ndim;

                if (otogen_is_silence(ctx->phdict, name)
                    || otogen_is_consonant(ctx->phdict, name))
                    continue;
                ds = dur_stats[stt];
                ds->ocpsum += 1.0f;
                ds->valsum += model->durations[i]->mean;
                ds->sqrsum += model->durations[i]->var;
                gs = gmm_stats[s];
                gs->ocpsum += 1.0f;
                gs->weightsum[0] += 1.0f;
                ndim = gs->ndim;
                if (ndim <= 0)
                    continue;
                gmm = stm->gmms[i];
                for (int k = 0; k < ndim; k++) {
                    gs->valsum[k] += gmm->mean[k];
                    gs->sqrsum[k] += gmm->var[k];
                }
            }
        }
    }

    for (int s = 0; s < nstream; s++) {
        lrh_gmm_stat *gs = gmm_stats[s];
        int ndim = gs->ndim;
        if (ndim <= 0)
            continue;
        for (int k = 0; k < ndim; k++) {
            gs->valsum[k] /= gs->ocpsum;
            gs->sqrsum[k] /= gs->ocpsum;
        }
    }
    for (int i = 0; i < nstate; i++) {
        lrh_duration_stat *ds = dur_stats[i];
        ds->valsum /= ds->ocpsum;
        ds->sqrsum = 1.0f;
    }

    if (ngmm > 0) {
        for (int s = 0; s < nstream; s++) {
            lrh_stream *stm = model->streams[s];
            lrh_gmm_stat *gs = gmm_stats[s];
            int ndim = gs->ndim;
            for (int i = 0; i < ngmm; i++) {
                int ph = i / nstate;
                wchar_t *name = ctx->phdict->names[ph];
                lrh_gmm *gmm;

                if (otogen_is_silence(ctx->phdict, name)
                    || otogen_is_consonant(ctx->phdict, name))
                    continue;
                if (ndim <= 0)
                    continue;
                gmm = stm->gmms[i];
                for (int k = 0; k < ndim; k++) {
                    gmm->mean[k] = gs->valsum[k];
                    gmm->var[k] = gs->sqrsum[k];
                }
            }
        }
    }

    for (int s = 0; s < nstream; s++)
        lrh_delete_gmm_stat(gmm_stats[s]);
    for (int i = 0; i < nstate; i++)
        lrh_delete_duration_stat(dur_stats[i]);
    free(gmm_stats);
    free(dur_stats);
}

void otogen_train(otogen_ctx *ctx, int niter, int hsmm)
{
    FP_TYPE prev = 0.0f;

    if (ctx->observset == NULL) {
        fwrite("Error: audio data not initialized yet.\n", 1, 0x27, stderr);
        return;
    }
    if (ctx->segset == NULL) {
        fwrite("Error: segmentation not initialized yet.\n", 1, 0x29, stderr);
        return;
    }
    if (ctx->model == NULL) {
        fwrite("Error: model not initialized yet.\n", 1, 0x22, stderr);
        return;
    }
    if (ctx->infer_nt == NULL) {
        fwrite("Error: list of inference settings not initialized yet.\n", 1, 0x37, stderr);
        return;
    }

    for (int it = 0; it < niter; it++) {
        lrh_model_stat *st;
        FP_TYPE sum = 0.0f;
        FP_TYPE mean;

        printf("Iteration %d/%d: ", it + 1, niter);
        fflush(stdout);
        st = lrh_model_stat_from_model(ctx->model);
        for (int i = 0; i < ctx->nfile; i++) {
            lrh_observ *ob = ctx->observset->samples[i];
            lrh_seg *seg = ctx->segset->samples[i];
            FP_TYPE ll;

            if (seg->nseg > 0)
                seg->time[seg->nseg - 1] = ob->nt;
            lrh_inference_stprune = (FP_TYPE)ctx->infer_nt[i];
            if (hsmm)
                ll = lrh_estimate(st, ctx->model, ob, seg);
            else
                ll = lrh_estimate_geometric(st, ctx->model, ob, seg);
            if (ll < -1000000000.0)
                fprintf(stderr, "Warning: unsuccessful inference on sample %d/%d.\n",
                    i + 1, ctx->nfile);
            sum += ll;
        }
        mean = sum / (FP_TYPE)ctx->nfile;
        printf("Mean log likelihood = %f\n", (double)mean);
        lrh_model_update(ctx->model, st, hsmm ? 0 : 1);
        lrh_delete_model_stat(st);
        if (it > 1 && prev + 1.0 >= (double)mean) {
            printf("Convergence reached at iteration %d.\n", it + 1);
            return;
        }
        prev = mean;
    }
}

void otogen_edit_silence_floor(otogen_ctx *ctx)
{

    if (ctx->model == NULL) {
        fwrite("Error: model not initialized yet.\n", 1, 0x22, stderr);
        return;
    }
    if (ctx->phdict == NULL) {
        fwrite("Error: phoneme dictionary not initialized yet.\n", 1, 0x2F, stderr);
        return;
    }
    for (int i = 0; i < ctx->phdict->all->count; i++) {
        wchar_t *name = ctx->phdict->names[i];
        if (!otogen_is_silence(ctx->phdict, name)
            && otogen_is_consonant(ctx->phdict, name)
            && ctx->nstate > 0) {
            for (int st = 0; st < ctx->nstate; st++)
                ctx->model->durations[i * ctx->nstate + st]->floor = ctx->nskip;
        }
    }
}

void otogen_align(otogen_ctx *ctx, int geometric)
{

    if (ctx->observset == NULL) {
        fwrite("Error: audio data not initialized yet.\n", 1, 0x27, stderr);
        return;
    }
    if (ctx->segset == NULL) {
        fwrite("Error: segmentation not initialized yet.\n", 1, 0x29, stderr);
        return;
    }
    if (ctx->model == NULL) {
        fwrite("Error: model not initialized yet.\n", 1, 0x22, stderr);
        return;
    }
    if (ctx->phdict == NULL) {
        fwrite("Error: phoneme dictionary not initialized yet.\n", 1, 0x2F, stderr);
        return;
    }
    if (ctx->infer_nt == NULL) {
        fwrite("Error: list of inference settings not initialized yet.\n", 1, 0x37, stderr);
        return;
    }

    for (int i = 0; i < ctx->nfile; i++) {
        lrh_observ *ob = ctx->observset->samples[i];
        lrh_seg *seg = ctx->segset->samples[i];
        FP_TYPE *outp;
        int *path;

        lrh_inference_stprune = (FP_TYPE)ctx->infer_nt[i];
        if (geometric) {
            outp = lrh_sample_outputprob(ctx->model, ob, seg);
            path = lrh_viterbi_geometric(ctx->model, seg, outp, ob->nt, NULL);
        } else {
            outp = lrh_sample_outputprob_lg(ctx->model, ob, seg);
            path = lrh_viterbi(ctx->model, seg, outp, ob->nt, NULL);
        }
        if (seg->nseg > 0) {
            for (int j = 0; j < seg->nseg; j++)
                seg->time[j] = path[j];
        }
        free(outp);
        free(path);

        if (ctx->labels != NULL && seg->nseg > 0) {
            otogen_label *lab = ctx->labels[i];
            FP_TYPE hop = ctx->hop_sec;
            int st = ctx->nstate;
            int p = 0;
            FP_TYPE prev = 0.0f;
            int idx = st;
            while (1) {
                lab->tbegin[p] = prev * hop;
                lab->tend[p] = hop * (FP_TYPE)seg->time[idx - 1];
                if (seg->nseg <= idx)
                    break;
                prev = (FP_TYPE)seg->time[idx - 1];
                idx += st;
                ++p;
            }
        }
    }
}

#include "otogen_boundary_model.inc"

otogen_bmodel *otogen_bmodel_from_table(void)
{
    const float *a1 = otogen_boundary_model_W;
    int nclass = 5;
    int ndim = 271;
    otogen_bmodel *m;
    int off;

    m = (otogen_bmodel *)malloc(sizeof(otogen_bmodel));
    m->nclass = nclass;
    m->ndim = ndim;
    m->cell = (otogen_bmodel_cell ***)malloc((size_t)nclass * sizeof(otogen_bmodel_cell **));
    for (int i = 0; i < nclass; i++) {
        m->cell[i] = (otogen_bmodel_cell **)malloc((size_t)nclass * sizeof(otogen_bmodel_cell *));
        for (int j = 0; j < nclass; j++) {
            otogen_bmodel_cell *c = (otogen_bmodel_cell *)malloc(sizeof(otogen_bmodel_cell));
            c->mean = (float *)calloc((size_t)ndim, sizeof(float));
            c->ndim = ndim;
            c->weight = 1e10f;
            m->cell[i][j] = c;
        }
    }
    off = 0;
    for (int i = 0; i < nclass; i++) {
        for (int j = 0; j < nclass; j++)
            m->cell[i][j]->weight = a1[off++];
    }
    {
        int mi = 0;
        for (int i = 0; i < nclass; i++) {
            for (int j = 0; j < nclass; j++) {
                if (m->cell[i][j]->weight >= 1.0f)
                    continue;
                for (int k = 0; k < ndim; k++)
                    m->cell[i][j]->mean[k] = otogen_boundary_model_mean[mi][k];
                mi++;
            }
        }
    }
    return m;
}

void otogen_bmodel_free(otogen_bmodel *m)
{

    if (m == NULL)
        return;
    for (int i = 0; i < m->nclass; i++) {
        for (int j = 0; j < m->nclass; j++) {
            free(m->cell[i][j]->mean);
            free(m->cell[i][j]);
        }
        free(m->cell[i]);
    }
    free(m->cell);
    free(m);
}

static int otogen_phoneme_class(otogen_phdict *d, const wchar_t *ph)
{
    if (otogen_is_voiced(d, ph))
        return 0;
    if (otogen_is_unvoiced(d, ph))
        return 1;
    if (otogen_is_semi(d, ph))
        return 2;
    if (otogen_is_vowel(d, ph))
        return 3;
    if (otogen_is_silence(d, ph))
        return 4;
    return -1;
}

static int otogen_fill_bfeat(FP_TYPE **band, int t, int nt, int nband, FP_TYPE *feat)
{
    FP_TYPE cep[8][10];

    if (t <= 3 || t + 9 >= nt)
        return 0;
    for (int i = 0; i < 4; i++) {
        FP_TYPE *left = cig_dct(band[t - 4 + i], nband);
        FP_TYPE *right = cig_dct(band[t + 1 + i], nband);
        for (int k = 0; k < 10; k++) {
            feat[i * 10 + k] = left[k];
            feat[40 + i * 10 + k] = right[k];
            cep[i][k] = left[k];
            cep[4 + i][k] = right[k];
        }
        free(left);
        free(right);
    }
    for (int i = 1; i < 8; i++) {
        for (int k = 0; k < 10; k++)
            feat[80 + (i - 1) * 10 + k] = (FP_TYPE)fabs((double)(cep[i][k] - cep[i - 1][k]));
    }
    for (int i = 0; i < 6; i++) {
        for (int k = 0; k < 10; k++) {
            FP_TYPE a = cep[i][k] + cep[i + 2][k] - 2.0f * cep[i + 1][k];
            feat[150 + i * 10 + k] = a * a;
        }
    }
    for (int i = 0; i < 6; i++) {
        for (int k = 0; k < 10; k++) {
            FP_TYPE a = cep[i + 2][k] - cep[i][k];
            feat[210 + i * 10 + k] = a * a;
        }
    }
    return 1;
}

static FP_TYPE otogen_lda_sigmoid(otogen_bmodel *bm, int pair, const FP_TYPE *feat)
{
    int c0 = pair / 5;
    int c1 = pair % 5;
    otogen_bmodel_cell *cell;
    const float *mean;
    FP_TYPE s;
    FP_TYPE p;

    if (bm == NULL)
        return 0.0000099999997f;
    cell = bm->cell[c0][c1];
    if (cell->weight >= 0.5f)
        return 0.0000099999997f;
    mean = cell->mean;
    s = mean[0];
    for (int k = 0; k < 270; k++)
        s += mean[k + 1] * feat[k];
    if (s < -20.0f)
        s = -20.0f;
    if (s > 20.0f)
        s = 20.0f;
    p = 1.0f / (exp_2((FP_TYPE)(-s)) + 1.0f);
    if (p <= 0.5f)
        p = 0.0000099999997f;
    return p;
}

void otogen_refine_alignment(otogen_ctx *ctx, int mode)
{
    int nband;

    if (ctx->observset == NULL) {
        fwrite("Error: audio data not initialized yet.\n", 1, 0x27, stderr);
        return;
    }
    if (ctx->phdict == NULL) {
        fwrite("Error: phoneme dictionary not initialized yet.\n", 1, 0x2F, stderr);
        return;
    }
    if (ctx->labels == NULL) {
        fwrite("Error: list of labels not initialized yet.\n", 1, 0x2B, stderr);
        return;
    }
    if (ctx->boundary_model == NULL) {
        fwrite("Error: model not initialized yet.\n", 1, 0x22, stderr);
        return;
    }

    nband = ctx->nmel_delta;
    for (int fi = 0; fi < ctx->nfile; fi++) {
        otogen_label *lab = ctx->labels[fi];
        lrh_observ *ob;
        FP_TYPE **band;
        int nt;
        FP_TYPE hop;
        FP_TYPE inv_hop;
        int win_big;
        int win_small;
        FP_TYPE **feat_at;
        FP_TYPE **prob_at;

        if (lab == NULL || ctx->band_e == NULL || ctx->band_e[fi] == NULL)
            continue;
        ob = ctx->observset->samples[fi];
        if (ob == NULL)
            continue;
        nt = ob->nt;
        band = ctx->band_e[fi];
        hop = ctx->hop_sec;
        inv_hop = 1.0f / hop;
        feat_at = (FP_TYPE **)calloc((size_t)nt, sizeof(FP_TYPE *));
        prob_at = (FP_TYPE **)calloc((size_t)nt, sizeof(FP_TYPE *));
        for (int t = 0; t < nt; t++) {
            if (t <= 3 || t + 9 >= nt)
                continue;
            feat_at[t] = (FP_TYPE *)calloc(270, sizeof(FP_TYPE));
            if (!otogen_fill_bfeat(band, t, nt, nband, feat_at[t])) {
                free(feat_at[t]);
                feat_at[t] = NULL;
                continue;
            }
            if (mode == 1) {
                prob_at[t] = (FP_TYPE *)calloc(25, sizeof(FP_TYPE));
                for (int pair = 0; pair < 25; pair++)
                    prob_at[t][pair] = otogen_lda_sigmoid(ctx->boundary_model, pair, feat_at[t]);
            }
        }

        if (mode == 1) {
            int nph = lab->nseg;
            int nseg = 2 * nph - 1;
            lrh_model *bm;
            lrh_seg *bseg;
            FP_TYPE *outp;
            int *path;

            bm = lrh_create_empty_model(0, nph + 1);
            for (int p = 0; p < nph; p++) {
                lrh_duration *d = lrh_create_duration();
                d->mean = (lab->tend[p] - lab->tbegin[p]) * inv_hop;
                d->var = 500.0f;
                d->floor = 3;
                bm->durations[p + 1] = d;
            }
            bm->durations[0] = lrh_create_duration();
            bm->durations[0]->mean = 0.0f;
            bm->durations[0]->var = 2.0f;
            bm->durations[0]->floor = 8;

            bseg = lrh_create_seg(0, nseg);
            for (int p = 0; p < nph; p++) {
                bseg->time[2 * p] = (int)(lab->tend[p] * inv_hop);
                bseg->durstate[2 * p] = p + 1;
                if (p != nph - 1) {
                    bseg->time[2 * p + 1] = (int)(lab->tend[p] * inv_hop + 1.0f);
                    bseg->durstate[2 * p + 1] = 0;
                }
            }

            outp = (FP_TYPE *)calloc((size_t)nt * (size_t)nseg, sizeof(FP_TYPE));
            {
                FP_TYPE floor_ll = log_2(ctx->silence_floor);
                for (int t = 0; t < nt; t++) {
                    for (int j = 0; j < nph; j++)
                        outp[t * nseg + 2 * j] = floor_ll;
                }
            }
            win_big = (int)(ctx->pad_sec * inv_hop);
            win_small = (int)(0.05f * inv_hop);
            for (int p = 0; p < nph - 1; p++) {
                int c0 = otogen_phoneme_class(ctx->phdict, lab->phoneme[p]);
                int c1 = otogen_phoneme_class(ctx->phdict, lab->phoneme[p + 1]);
                int tb = (int)round((double)(lab->tbegin[p] * inv_hop));
                int tm = (int)round((double)(lab->tend[p] * inv_hop));
                int te = (int)round((double)(lab->tend[p + 1] * inv_hop));
                int lo;
                int hi;
                int pair = (c0 >= 0 && c1 >= 0) ? c0 * 5 + c1 : -1;
                for (int t = 0; t < nt; t++)
                    outp[t * nseg + 2 * p + 1] = -1.0e10f;
                lo = tb - win_small;
                if (lo < tm - win_big)
                    lo = tm - win_big;
                if (lo < 0)
                    lo = 0;
                hi = te + win_small;
                if (hi > tm + win_big)
                    hi = tm + win_big;
                if (hi > nt)
                    hi = nt;
                if (pair >= 0) {
                    for (int t = lo; t < hi; t++) {
                        if (prob_at[t] != NULL)
                            outp[t * nseg + 2 * p + 1] = prob_at[t][pair];
                    }
                }
            }
            path = lrh_viterbi(bm, bseg, outp, nt, NULL);
            if (path != NULL && nph > 1) {
                for (int p = 1; p < nph; p++) {
                    FP_TYPE mid = 0.5f * hop * (FP_TYPE)(path[2 * p - 2] + path[2 * p - 1]);
                    lab->tend[p - 1] = mid;
                    lab->tbegin[p] = mid;
                }
            }
            free(path);
            lrh_delete_seg(bseg);
            lrh_delete_model(bm);
            free(outp);
            for (int i = 0; i < nt; i++)
                free(prob_at[i]);
        } else if (lab->nseg > 2) {
            win_big = (int)(0.2f * inv_hop);
            win_small = (int)(0.05f * inv_hop);
            for (int p = 1; p < lab->nseg - 1; p++) {
                int c0 = otogen_phoneme_class(ctx->phdict, lab->phoneme[p - 1]);
                int c1 = otogen_phoneme_class(ctx->phdict, lab->phoneme[p]);
                int tb;
                int tm;
                int te;
                int lo;
                int hi;
                FP_TYPE best;
                int best_t;

                if (c0 < 0 || c1 < 0)
                    continue;
                int pair = c0 * 5 + c1;
                if (ctx->boundary_model->cell[c0][c1]->weight > 0.5f)
                    continue;
                if ((c0 | c1) == 0
                    && otogen_is_consonant(ctx->phdict, lab->phoneme[p - 1])
                    && otogen_is_consonant(ctx->phdict, lab->phoneme[p]))
                    continue;
                tb = (int)round((double)(lab->tbegin[p - 1] * inv_hop));
                tm = (int)round((double)(lab->tbegin[p] * inv_hop));
                te = (int)round((double)(lab->tend[p] * inv_hop));
                lo = tb - win_small;
                if (lo < tm - win_big)
                    lo = tm - win_big;
                if (lo < 0)
                    lo = 0;
                hi = te + win_small;
                if (hi > tm + win_big)
                    hi = tm + win_big;
                if (hi > nt)
                    hi = nt;
                best = 0.0f;
                best_t = tm;
                for (int t = lo; t < hi; t++) {
                    FP_TYPE s;
                    if (feat_at[t] == NULL)
                        continue;
                    s = otogen_lda_sigmoid(ctx->boundary_model, pair, feat_at[t]);
                    if (s > best) {
                        best = s;
                        best_t = t;
                    }
                }
                if (best <= 0.5f)
                    best_t = tm;
                lab->tend[p - 1] = hop * (FP_TYPE)best_t;
                lab->tbegin[p] = hop * (FP_TYPE)best_t;
            }
        }

        for (int t = 0; t < nt; t++)
            free(feat_at[t]);
        free(feat_at);
        free(prob_at);
    }
}

void otogen_fill_label_times(otogen_ctx *ctx)
{

    if (ctx->nfile <= 0)
        return;
    for (int i = 0; i < ctx->nfile; i++) {
        otogen_label *lab = ctx->labels[i];
        FP_TYPE pad = ctx->pad_sec;
        FP_TYPE t0 = ctx->trim[2 * i];
        FP_TYPE t1 = ctx->trim[2 * i + 1];
        FP_TYPE dur;

        dur = (FP_TYPE)ctx->observset->samples[i]->nt * ctx->hop_sec - (pad + pad);
        for (int k = 0; k < lab->nseg; k++) {
            FP_TYPE a = lab->tbegin[k] + t0 - pad;
            FP_TYPE b = lab->tend[k] - pad;
            if (a < 0.0f)
                a = 0.0f;
            if (dur < b)
                b = dur;
            lab->tbegin[k] = a;
            lab->tend[k] = b + t0;
        }
        lab->tend[lab->nseg - 1] += t1;
        lab->tbegin[0] = 0.0f;
    }
}
