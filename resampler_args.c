#include "resampler.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#include "external/mrq/mrq.h"
#include "flag.h"
#include "logging.h"


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
        return (unsigned char)(c - 'A');
    if ((unsigned char)(c - 'a') <= 25u)
        return (unsigned char)(c - 71);
    if ((unsigned char)(c - '0') <= 9u)
        return (unsigned char)(c + 4);
    return (c != '+') + 62;
}

// Decode the upper sextet in its final bit position.
static int pb_base64_high(unsigned char c)
{
    if ((unsigned char)(c - 'A') <= 25u)
        return (unsigned char)(c - 'A') * 64;
    if ((unsigned char)(c - 'a') <= 25u)
        return (unsigned char)(c - 71) * 64;
    if ((unsigned char)(c - '0') <= 9u)
        return (unsigned char)(c + 4) * 64;
    return c == '+' ? 62 * 64 : 63 * 64;
}

static void parse_pitchbend(resampler_args *args, const char *src)
{
    char *pb = strdup(src);
    args->pitch_bend = calloc(20, sizeof(float));
    int len = (int)strlen(pb);
    int cap = 20;
    int i = 0;

    for (; i < len; ++i) {
        if (pb[i] == '#') {
            int start = i;
            while (i < len - 1) {
                if (pb[++i] == '#')
                    break;
            }
            pb[i] = 0;
            int count = atoi(&pb[start + 1]);
            int old_size = args->pitch_bend_size;
            int new_size = old_size + count;
            if (cap < new_size) {
                cap += count;
                args->pitch_bend = realloc(args->pitch_bend,
                                           (size_t)cap * sizeof(float));
                old_size = args->pitch_bend_size;
            }
            if (count > 0) {
                float *previous = args->pitch_bend + old_size - 1;
                float *cursor = args->pitch_bend + old_size;
                float *limit = args->pitch_bend + new_size;
                do {
                    *cursor++ = *previous;
                } while (cursor != limit);
            }
            args->pitch_bend_size = new_size;
        } else {
            int value = pb_base64_high((unsigned char)pb[i]) +
                        pb_base64_value((unsigned char)pb[i + 1]);
            if (value >= 2048)
                value -= 4096;
            if (cap <= args->pitch_bend_size) {
                cap = args->pitch_bend_size * 2;
                args->pitch_bend = realloc(args->pitch_bend,
                                           (size_t)cap * sizeof(float));
            }
            ++i;
            args->pitch_bend[args->pitch_bend_size++] = (float)value;
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

    if (argc > 5) {
        args->flags = strdup(argv[5]);

        flag_table *f = malloc(sizeof(*f));
        f->capacity = FLG_COUNT;
        f->values = calloc(FLG_COUNT, sizeof(float));
        f->set_flags = calloc(FLG_COUNT, sizeof(int));
        f->values[FLG_Mm] = 100.0f;
        f->values[FLG_P] = 86.0f;
        flag_parse(f, args->flags);

        float *v = f->values;
        int *set = f->set_flags;

        args->g = clampf(v[FLG_g], -100.0f, 100.0f);
        args->t = clampf(v[FLG_t], -1200.0f, 1200.0f);
        args->Y = v[FLG_Y];
        if (set[FLG_e] > 0)
            args->e = 1;
        else
            args->e = -(set[FLG_Me] > 0);
        args->b = fmax(-19.9999999999, fminf(v[FLG_b], 100.0f));
        args->A = clampf(v[FLG_A], -100.0f, 100.0f);
        args->P = clampf(v[FLG_P], 0.0f, 100.0f);
        args->u = set[FLG_u] > 0;
        args->Mt = clampf(v[FLG_Mt], -100.0f, 100.0f);
        args->Mb = clampf(v[FLG_Mb], -100.0f, 100.0f);
        args->Md = clampf(v[FLG_Md], -100.0f, 100.0f);
        args->Mm = clampf(v[FLG_Mm], 0.0f, 100.0f);
        {
            double ms = v[FLG_Ms];
            double rounded = round(ms);
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
        args->_s = f->values[FLG__s];

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
                                }
                                if (argc > 13 && args->tempo > 0.0f) {
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
