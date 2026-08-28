#include "resampler.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#include "config.h"
#include "external/ciglet/ciglet.h"
#include "external/libllsm/llsm.h"
#include "external/mrq/mrq.h"
#include "logging.h"
#include "utils.h"
#include "flag.h"
#include "lmix.h"
#include "xchg.h"



static float clampf(float x, float lo, float hi)
{
    if (x > hi)
        x = hi;
    if (x < lo)
        x = lo;
    return x;
}

// Decode one UTAU pitchbend alphabet character.
static int pb_base64_value(unsigned char c)
{
    if ((unsigned char)(c - 'A') <= 25u)
        return c - 'A';
    if ((unsigned char)(c - 'a') <= 25u)
        return (unsigned char)(c - 71);
    if ((unsigned char)(c - '0') <= 9u)
        return (unsigned char)(c + 4);
    return (c != '+') + 62;
}

static void parse_pitchbend(resampler_args *args, const char *src)
{
    char *pb = strdup(src);
    int len = (int)strlen(pb);
    int cap = 20;
    int i = 0;

    args->pitch_bend = calloc(20, sizeof(float));
    while (i < len) {
        if (pb[i] == '#') {
            int start = i;
            while (i < len - 1) {
                if (pb[++i] == '#')
                    break;
            }
            pb[i] = 0;
            int count = atoi(&pb[start + 1]);
            int new_size = args->pitch_bend_size + count;
            if (cap < new_size) {
                cap += count;
                args->pitch_bend = realloc(args->pitch_bend,
                                           (size_t)cap * sizeof(float));
            }
            if (count > 0) {
                float last = args->pitch_bend[args->pitch_bend_size - 1];
                for (int k = 0; k < count; k++)
                    args->pitch_bend[args->pitch_bend_size + k] = last;
            }
            i++;
            args->pitch_bend_size = new_size;
        } else {
            int value = pb_base64_value((unsigned char)pb[i]) * 64 +
                        pb_base64_value((unsigned char)pb[i + 1]);
            if (value >= 2048)
                value -= 4096;
            if (cap <= args->pitch_bend_size) {
                cap = args->pitch_bend_size * 2;
                args->pitch_bend = realloc(args->pitch_bend,
                                           (size_t)cap * sizeof(float));
            }
            args->pitch_bend[args->pitch_bend_size++] = (float)value;
            i += 2;
        }
    }
    free(pb);
}

resampler_args *parse_resampler_args(int argc, char **argv)
{
    if (argc <= 4)
        return NULL;

    resampler_args *args = calloc(1, sizeof(resampler_args));
    args->in_path = malloc(0x800);
    args->out_path = malloc(0x800);
    args->error_code = 0;
    args->skip_llsm = 0;

    int wide_argc = argc;
    LPWSTR *wide_argv = CommandLineToArgvW(GetCommandLineW(), &wide_argc);
    args->in_path = wcscpy(args->in_path, wide_argv[1]);
    args->out_path = wcscpy(args->out_path, wide_argv[2]);
    if (wcscmp(args->out_path, L"nul") == 0)
        args->error_code = 1;

    args->pitch = strdup(argv[3]);
    args->velocity = (float)atof(argv[4]);
    args->length_require = -1.0f;
    args->volume = 100.0f;
    args->modulation = 0.0f;
    args->unknow12 = 0.0f;

    if (argc != 5) {
        args->flags = strdup(argv[5]);

        flag_table *f = flag_create();
        f->values[FLG_Mm] = 100.0f;
        f->values[FLG_P] = 86.0f;
        flag_parse(100.0f, f, args->flags);

        float *v = f->values;
        int *set = f->set_flags;

        args->g = clampf(v[FLG_g], -100.0f, 100.0f);
        args->t = clampf(v[FLG_t], -1200.0f, 1200.0f);
        args->Y = v[FLG_Y];
        if (set[FLG_e] > 0)
            args->e = 1;
        else
            args->e = -(set[FLG_Me] > 0);
        args->b = clampf(v[FLG_b], -20.0f, 100.0f);
        args->A = clampf(v[FLG_A], -100.0f, 100.0f);
        args->P = clampf(v[FLG_P], 0.0f, 100.0f);
        args->u = set[FLG_u] > 0;
        args->Mt = clampf(v[FLG_Mt], -100.0f, 100.0f);
        args->Mb = clampf(v[FLG_Mb], -100.0f, 100.0f);
        args->Md = clampf(v[FLG_Md], -100.0f, 100.0f);
        args->Mm = clampf(v[FLG_Mm], 0.0f, 100.0f);
        {
            double ms = v[FLG_Ms];
            double rounded = rint(fabs(ms) + 0.5);
            if (ms < 0.0)
                rounded = -rounded;
            if (rounded > 10.0)
                rounded = 10.0;
            if (rounded < 0.0)
                rounded = 0.0;
            args->Ms = (int)rounded;
        }
        args->ME = clampf(v[FLG_ME], -100.0f, 100.0f);
        args->Mo = clampf(v[FLG_Mo], -100.0f, 100.0f);
        args->Mr = clampf(v[FLG_Mr], -100.0f, 100.0f);
        args->MC = clampf(v[FLG_MC], 0.0f, 100.0f);
        args->MG = clampf(v[FLG_MG], 0.0f, 100.0f);
        args->MD = clampf(v[FLG_MD], 0.0f, 100.0f);
        args->MP = clampf(v[FLG_Mp], 0.0f, 100.0f);
        if (set[FLG__m] > 0)
            args->modulation = v[FLG__m];
        if (set[FLG__f] > 0)
            args->skip_llsm = 1;
        if (set[FLG__o] > 0)
            mrq_set_version((int)v[FLG__o]);
        args->_s = v[FLG__s];

        free(f->values);
        free(f->set_flags);
        free(f);

        if (argc > 7) {
            args->offset = (float)(atof(argv[6]) * 0.001);
            args->length_require = (float)(atof(argv[7]) * 0.001);
            if (argc != 8) {
                args->fixed_length = (float)(atof(argv[8]) * 0.001);
                if (argc != 9) {
                    args->end_blank = (float)(atof(argv[9]) * 0.001);
                    if (argc != 10) {
                        args->volume = (float)atof(argv[10]);
                        if (argc != 11) {
                            args->modulation = (float)atof(argv[11]);
                            if (argc != 12) {
                                args->tempo = (float)atof(argv[12] + 1);
                                if (args->tempo <= 0.0f) {
                                    log_file_printf(
                                        "Warning: invalid BPM sequence '%ls'.\n",
                                        (const wchar_t *)argv[12]);
                                } else if (argc != 13) {
                                    parse_pitchbend(args, argv[13]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (args->pitch_bend_size > 0) {
        args->pitch_bend_size++;
        args->pitch_bend = realloc(args->pitch_bend,
                                   (size_t)args->pitch_bend_size * sizeof(float));
        args->pitch_bend[args->pitch_bend_size - 1] = 0.0f;
    } else {
        args->pitch_bend = calloc(1, sizeof(float));
        args->pitch_bend_size = 1;
        args->pitch_bend[0] = 0.0f;
    }
    return args;
}

void free_resampler_args(resampler_args *args)
{
    if (args != NULL) {
        if (args->in_path)
            free(args->in_path);
        if (args->out_path)
            free(args->out_path);
        if (args->pitch)
            free(args->pitch);
        if (args->flags)
            free(args->flags);
        if (args->pitch_bend)
            free(args->pitch_bend);
        free(args);
    }
}

int main_resampler(int argc, char **argv)
{
    resampler_args *args = parse_resampler_args(argc, argv);
    log_file_printf("main_resampler: input & output: %s, %s\n",
                    argv[1], argv[2]);
    log_file_printf("pitch percent: %s, velocity: %f\n",
                    args->pitch, args->velocity);
    log_file_printf("flags: %s\n", args->flags);
    log_file_printf(
        "offset: %f, length require: %f, fixed length: %f, end blank: %f\n",
        args->offset, args->length_require, args->fixed_length, args->end_blank);
    log_file_printf("volume: %f\n", args->volume);
    log_file_printf("modulation: %f, bpm: %f, size of pitch bend: %d\n",
                    args->modulation, args->tempo, args->pitch_bend_size);

    wchar_t *in_dir = _wcsdup(args->in_path);
    *path_get_filename(in_dir) = 0;
    wchar_t *cfg_path = wstr_concat(in_dir, L"moreconfig.txt");
    load_config(cfg_path, 0);
    free(in_dir);
    free(cfg_path);

    lmix_data *lmix;
    if (args->u != 0) {
        lmix = lmix_from_file_direct(args->in_path);
    } else {
        omp_set_num_threads(omp_get_max_threads());
        float end_len;
        if (args->end_blank < 0.0f)
            end_len = args->offset - args->end_blank;
        else
            end_len = -args->end_blank;
        float start = args->offset - 0.05f;
        if (start < 0.0f)
            start = 0.0f;
        lmix = lmix_from_file(args->in_path, start, end_len, args->skip_llsm);
        if (g_multithread_synthesis == 0)
            omp_set_num_threads(1);
    }

    if (lmix == NULL) {
        log_file_printf(
            "main_resampler: cannot open llsm file. Synthesis terminated.\n");
        return 1;
    }

    if (args->error_code != 0) {
        lmix_delete(lmix);
        free_resampler_args(args);
        log_file_printf(
            "main_resampler: analysis-only resampler session completed.\n");
        return 0;
    }

    float total_dur = (float)lmix->conf.nfrm * lmix->conf.thop;
    float stretch_end;
    if (args->end_blank < 0.0f)
        stretch_end = args->end_blank - args->offset + total_dur;
    else
        stretch_end = args->end_blank;
    float stretch = args->fixed_length *
                    (float)pow(2.0, 1.0 - 0.01 * args->velocity);
    float length = total_dur;
    if (args->length_require >= 0.0f)
        length = args->length_require;

    log_file_printf("main_resampler: modifying parameters...\n");
    lmix_data *lmix2;
    if (args->u != 0) {
        lmix2 = lmix_create(lmix->conf);
        for (int i = 0; i < lmix2->conf.nfrm; i++) {
            lmix2->frames[i] = lmix_create_frame(lmix->conf, g_lmix_nfft, 0);
            lmix_copy_frame(lmix2->frames[i], lmix->frames[i]);
        }
        llsm_wave *wave = lmix2->frames[0]->wave;
        if (wave != NULL) {
            int skip = (int)rint((double)wave->fs * (args->offset - 0.01));
            if (skip > wave->nsample)
                skip = wave->nsample;
            if (skip < 0)
                skip = 0;
            int remain = wave->nsample - skip;
            if (remain > 0) {
                for (int j = 0; j < remain; j++)
                    wave->samples[j] = wave->samples[skip + j];
            }
            wave->nsample = remain;
            wave->t0 = 0.01f;
        }
    } else {
        lmix2 = lmix_modify_data(lmix, args, args->offset, stretch_end,
                                 args->fixed_length, total_dur, stretch, length,
                                 0, NULL, NULL);
    }

    if (g_synthesis_utau_style_normalization != 0) {
        log_file_printf("main_resampler: normalizing...\n");
        lmix_normalize_data(lmix2, 0.005f * args->volume, args->P);
    } else {
        lmix_scale_volume(lmix2, 0.01f * args->volume);
    }

    if (args->Mb > 0.0f || args->MC > 0.0f) {
        log_file_printf("main_resampler: applying 'Mb' and/or 'MC' flag...\n");
        float mbmc = args->MC;
        if (mbmc < args->Mb)
            mbmc = args->Mb;
        float gain = (float)log(1.0000000001 - mbmc * 0.01);
        for (int i = 0; i < lmix2->conf.nfrm; i++) {
            lmix_frame *frame = lmix2->frames[i];
            int n = frame->nhar * frame->vt_stride;
            for (int k = 0; k < n; k++)
                frame->vt_magn[k] += gain;
        }
    }

    if (g_resampler_compatibility != 0) {
        int fs = g_output_sampling_rate;

        llsm_parameters param = llsm_init(lmix2->conf.nnosband);
        log_file_printf("main_resampler: converting LMIX to lv0 LLSM...\n");
        llsm_layer0 *layer0 = lmix_to_layer0(lmix2);
        log_file_printf("main_resampler: synthesizing...\n");
        param.s_fs = fs;
        llsm_output *out = llsm_layer0_synthesize(param, layer0);
        FP_TYPE *wave = lmix_render_wave(out, lmix2);
        llsm_deinit(param);
        llsm_delete_layer0(layer0);
        log_file_printf("main_resampler: outputing wav file...\n");
        char *out_mbs = wchar_to_mbs(args->out_path);
        wavwrite(wave, out->ny, fs, g_output_bit_depth, out_mbs);
        free(out_mbs);
        llsm_delete_output(out);
        free(wave);
    } else {
        log_file_printf("main_resampler: creating dummy wav file...\n");
        FILE *fp = _wfopen(args->out_path, L"w");
        fwrite("This is a dummy file created just to fool UTAU's cache mechanism "
               "and you really shouldn't expect to see anything here. The vocals "
               "are stored in llsm files, if nothing went wrong.\n",
               1, 0xB3u, fp);
        fclose(fp);
    }

    {
        log_tee_printf("main_resampler: outputing LMIX data pack...\n");
        wchar_t *tmp_path = calloc(wcslen(args->out_path) + 10, sizeof(wchar_t));
        wcscpy(tmp_path, args->out_path);
        wcscat(tmp_path, L".llsm.tmp");
        xchg_node *node = lmix_to_xchg(lmix2);
        log_file_printf("output_lmix: xchg_node created.\n");
        xchg_export_counted(node, tmp_path);
        log_file_printf("output_lmix: data file exported.\n");
        xchg_node_destroy(node);
        free(tmp_path);
    }

    log_file_printf("main_resampler: freeing memory...\n");
    lmix_delete(lmix);
    lmix_delete(lmix2);
    free_resampler_args(args);
    log_file_printf("main_resampler: resampler session completed.\n");
    return 0;
}
