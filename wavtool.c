#include "wavtool.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#include "config.h"
#include "external/ciglet/ciglet.h"
#include "external/libllsm/llsm.h"
#include "logging.h"
#include "utils.h"
#include "lmix.h"
#include "xchg.h"

static int lmix_conf_compatible(const llsm_conf *a, const llsm_conf *b)
{
    if (a->nhop != b->nhop) {
        log_tee_printf(
            "Non-conformant hop size between models to be concatenated.\n");
        return 0;
    }
    if (a->nnos != b->nnos) {
        log_tee_printf(
            "Non-conformant spectrum size between models to be concatenated.\n");
        return 0;
    }
    if (a->noswarp != b->noswarp) {
        log_tee_printf(
            "Non-conformant warping curve between models to be concatenated.\n");
        return 0;
    }
    if (a->nosf != b->nosf) {
        log_tee_printf(
            "Non-conformant frequency upperbound between models to be concatenated.\n");
        return 0;
    }
    if (a->nnosband != b->nnosband) {
        log_tee_printf(
            "Non-conformant band order between models to be concatenated.\n");
        return 0;
    }
    return 1;
}

// Convert wavtool arguments to an lmix envelope.
static lmix_envelope wavtool_envelope_from_args(const wavtool_args *args)
{
    // Moresampler uses (const lmix_envelope *)args to convert wavtool parameters to lmix envelopes.
    // This works in a 32-bit environment, but causes an alignment issue in a 64-bit environment;
    // therefore, a safe conversion is used here.
    lmix_envelope env = {0};

    env.t0 = args->offset;
    env.duration = args->length;
    env.p1 = args->p1;
    env.p2 = args->p2;
    env.p3 = args->p3;
    env.p4 = args->p4;
    env.p5 = args->p5;
    env.v1 = args->v1;
    env.v2 = args->v2;
    env.v3 = args->v3;
    env.v4 = args->v4;
    env.v5 = args->v5;
    return env;
}

static lmix_data *wavtool_empty_lmix(const llsm_conf *src, int nfrm)
{
    llsm_conf conf = *src;
    conf.nfrm = nfrm;
    lmix_data *out = lmix_create(conf);
    for (int i = 0; i < nfrm; ++i)
        out->frames[i] = lmix_create_frame(conf, g_lmix_nfft, 0);
    return out;
}

void free_wavtool_args(wavtool_args *args)
{
    if (args == NULL)
        return;
    if (args->input_path != NULL)
        free(args->input_path);
    if (args->output_path != NULL)
        free(args->output_path);
    free(args);
}

wavtool_args *parse_wavtool_args(int argc, char **argv)
{
    if (argc <= 4)
        return NULL;

    wavtool_args *args = (wavtool_args *)calloc(1, sizeof(wavtool_args));
    args->input_path = (wchar_t *)malloc(0x800);
    args->output_path = (wchar_t *)malloc(0x800);
    int narg = argc;
    LPWSTR *wargv = CommandLineToArgvW(GetCommandLineW(), &narg);
    wcscpy(args->output_path, wargv[1]);
    wcscpy(args->input_path, wargv[2]);
    args->last_note = 0;
    args->offset = (float)(atof(argv[3]) * 0.001);

    const char *len_str = argv[4];
    int len_n = (int)strlen(len_str);
    int has_at = 0;
    if (len_n > 0) {
        for (int i = 0; i < len_n; ++i) {
            if (len_str[i] == '@')
                has_at = 1;
        }
    }
    if (has_at) {
        char *dup = strdup(len_str);
        char *at = dup;
        int skip = 1;

        if (dup[0] == '@') {
            at = dup;
            skip = 1;
        } else {
            char *p = dup + 1;
            skip = 0;
            do {
                at = p++;
                ++skip;
            } while (p[-1] != '@');
            skip += 1;
        }
        *at = 0;
        char *rest = dup + skip;
        double head = atof(dup);
        char ch = *rest;
        char *signp = rest;
        int nsign = 1;
        if (((ch - '+') & 0xFD) != 0) {
            char *q = rest + 1;
            nsign = 0;
            do {
                ch = *q;
                ++nsign;
                signp = q++;
            } while (((ch - '+') & 0xFD) != 0);
            nsign += 1;
        } else {
            signp = rest;
            nsign = 1;
        }
        *signp = 0;
        double sign = (ch == '+') ? 1.0 : -1.0;
        double tempo = atof(rest);
        double corr = atof(rest + nsign);
        args->length = (float)(corr * sign * 0.001 + head * 0.125 / tempo);
        free(dup);
    } else {
        args->length = (float)(atof(len_str) * 0.001);
    }

    args->p1 = 0.0f;
    args->p2 = 0.0f;
    args->p3 = 0.0f;
    args->v1 = 100.0f;
    args->v2 = 100.0f;
    args->v3 = 100.0f;
    args->v4 = 100.0f;
    args->overlap = 0.0f;
    float v2 = 100.0f;

    if (argc == 5) {
        args->p4 = 0.0f;
        args->p5 = 0.0f;
        args->v5 = v2;
        if (strcmp(argv[argc - 1], "LAST_NOTE") == 0)
            args->last_note = 1;
        return args;
    }

    args->p1 = (float)(atof(argv[5]) * 0.001);
    if (argc == 6) {
        args->p4 = 0.0f;
        args->p5 = 0.0f;
        args->v5 = args->v2;
        if (strcmp(argv[5], "LAST_NOTE") == 0) {
            args->last_note = 1;
            args->p1 = 0.0f;
        }
        return args;
    }

    args->p2 = (float)(atof(argv[6]) * 0.001);
    if (argc == 7) {
        args->p4 = 0.0f;
        args->p5 = 0.0f;
        args->v5 = args->v2;
        if (strcmp(argv[6], "LAST_NOTE") == 0) {
            args->last_note = 1;
            args->p2 = 0.0f;
        }
        return args;
    }

    args->p3 = (float)(atof(argv[7]) * 0.001);
    if (argc == 8)
        goto fill_tail;

    args->v1 = (float)atof(argv[8]);
    if (argc == 9)
        goto fill_tail;

    v2 = (float)atof(argv[9]);
    args->v2 = v2;
    if (argc == 10)
        goto fill_tail;

    args->v3 = (float)atof(argv[10]);
    if (argc == 11)
        goto fill_tail;

    args->v4 = (float)atof(argv[11]);
    if (argc == 12)
        goto fill_tail;

    args->overlap = (float)(atof(argv[12]) * 0.001);
    args->p4 = 0.0f;
    args->p5 = 0.0f;
    args->v5 = args->v2;
    if (argc == 13) {
        if (strcmp(argv[argc - 1], "LAST_NOTE") == 0) {
            args->last_note = 1;
            args->overlap = 0.0f;
        }
        return args;
    }

    args->p4 = (float)(atof(argv[13]) * 0.001);
    if (argc != 14) {
        args->p5 = (float)(atof(argv[14]) * 0.001);
        if (argc != 15)
            args->v5 = (float)atof(argv[15]);
    }
    if (strcmp(argv[argc - 1], "LAST_NOTE") != 0)
        return args;
    args->last_note = 1;
    switch (argc) {
    case 9:
        args->v1 = 100.0f;
        return args;
    case 10:
        args->v2 = 100.0f;
        return args;
    case 11:
        args->v3 = 100.0f;
        return args;
    case 12:
        args->v4 = 100.0f;
        return args;
    case 13:
        args->overlap = 0.0f;
        return args;
    case 14:
        args->p4 = 0.0f;
        break;
    case 15:
        args->p5 = 0.0f;
        break;
    case 16:
        args->v5 = args->v2;
        break;
    default:
        break;
    }
    return args;

fill_tail:
    args->p4 = 0.0f;
    args->p5 = 0.0f;
    args->v5 = v2;
    if (strcmp(argv[argc - 1], "LAST_NOTE") != 0)
        return args;
    args->last_note = 1;
    if (argc == 8) {
        args->p3 = 0.0f;
        return args;
    }
    switch (argc) {
    case 9:
        args->v1 = 100.0f;
        return args;
    case 10:
        args->v2 = 100.0f;
        return args;
    case 11:
        args->v3 = 100.0f;
        return args;
    case 12:
        args->v4 = 100.0f;
        return args;
    case 13:
        args->overlap = 0.0f;
        return args;
    default:
        break;
    }
    return args;
}

static wavtool_args **wavtool_load_index(wchar_t *path, int *count)
{
    FILE *fp = _wfopen(path, L"rb");
    if (fp == NULL) {
        log_tee_printf("Error: cannot open indexing file %ls\n", path);
        return NULL;
    }
    xchg_file_lock(fp);
    wavtool_args **list = (wavtool_args **)calloc(10, sizeof(wavtool_args *));
    *count = 0;
    int cap = 10;
    wchar_t ch = 0;
    for (;;) {
        int n = 0;
        wchar_t *outp = (wchar_t *)calloc(0x400, sizeof(wchar_t));
        for (;;) {
            if (fread(&ch, 2, 1, fp) == 0) {
                free(outp);
                fclose(fp);
                return list;
            }
            if (ch == 0)
                break;
            outp[n++] = ch;
        }
        outp[n] = 0;
        if (outp == NULL)
            break;
        if (cap <= *count + 1) {
            cap *= 2;
            list = (wavtool_args **)realloc(list, (size_t)cap * sizeof(*list));
        }
        wavtool_args *ent = (wavtool_args *)malloc(sizeof(wavtool_args));
        ent->output_path = outp;
        list[*count] = ent;
        int m = 0;
        wchar_t *inp = (wchar_t *)calloc(0x400, sizeof(wchar_t));
        for (;;) {
            if (fread(&ch, 2, 1, fp) == 0) {
                free(inp);
                inp = NULL;
                break;
            }
            if (ch == 0)
                break;
            inp[m++] = ch;
        }
        if (inp != NULL)
            inp[m] = 0;
        ent->input_path = inp;
        fread(&ent->offset, 4, 1, fp);
        fread(&ent->length, 4, 1, fp);
        fread(&ent->overlap, 4, 1, fp);
        fread(&ent->p1, 4, 1, fp);
        fread(&ent->p2, 4, 1, fp);
        fread(&ent->p3, 4, 1, fp);
        fread(&ent->p4, 4, 1, fp);
        fread(&ent->p5, 4, 1, fp);
        fread(&ent->v1, 4, 1, fp);
        fread(&ent->v2, 4, 1, fp);
        fread(&ent->v3, 4, 1, fp);
        fread(&ent->v4, 4, 1, fp);
        fread(&ent->v5, 4, 1, fp);
        ++*count;
    }
    fclose(fp);
    return list;
}

static int wavtool_append_index(wchar_t *path, wavtool_args *args)
{
    FILE *fp = _wfopen(path, L"rb+");
    xchg_file_lock(fp);
    fseek(fp, 0, SEEK_END);
    fwrite(args->output_path, 2, wcslen(args->output_path), fp);
    wchar_t z = 0;
    fwrite(&z, 2, 1, fp);
    fwrite(args->input_path, 2, wcslen(args->input_path), fp);
    fwrite(&z, 2, 1, fp);
    fwrite(&args->offset, 4, 1, fp);
    fwrite(&args->length, 4, 1, fp);
    fwrite(&args->overlap, 4, 1, fp);
    fwrite(&args->p1, 4, 1, fp);
    fwrite(&args->p2, 4, 1, fp);
    fwrite(&args->p3, 4, 1, fp);
    fwrite(&args->p4, 4, 1, fp);
    fwrite(&args->p5, 4, 1, fp);
    fwrite(&args->v1, 4, 1, fp);
    fwrite(&args->v2, 4, 1, fp);
    fwrite(&args->v3, 4, 1, fp);
    fwrite(&args->v4, 4, 1, fp);
    fwrite(&args->v5, 4, 1, fp);
    return fclose(fp);
}

static int check_temp_last(wchar_t *bat_path, wchar_t *note_name)
{
    FILE *fp = _wfopen(bat_path, L"r");
    if (fp == NULL) {
        log_file_printf("check_temp_last: cannot open temp directory.\n");
        return 0;
    }
    fseek(fp, 0, SEEK_END);
    int size = ftell(fp);
    int window = size <= 1000 ? size : 1000;
    int off = size - window;
    if (size < window) {
        log_file_printf("check_temp_last: cannot locate the last note.\n");
        return 0;
    }
    wchar_t *last = NULL;
    do {
        int stop = off + 1200;

        fseek(fp, off, SEEK_SET);
        for (;;) {
            wchar_t *line = fgets_utf8_to_wchar(fp);

            if (line == NULL)
                break;
            if (wcsncmp(line, L"@set temp=", 10) == 0) {
                int n = (int)wcslen(line) - 1;
                int k = n;
                wchar_t *base = line;

                if (k < 0) {
                    base = line;
                } else {
                    wchar_t c = line[k];
                    if (c != L'/' ) {
                        while (c != L'\\') {
                            if (--k == -1) {
                                base = line;
                                goto got_base;
                            }
                            c = line[k];
                            n = k;
                            if (c == L'/')
                                break;
                        }
                    }
                    base = &line[n + 1];
                }
got_base:
                if (last != NULL)
                    free(last);
                last = _wcsdup(base);
                if (ftell(fp) > stop)
                    break;
            } else if (ftell(fp) > stop) {
                break;
            }
            free(line);
        }
        window += 1000;
        off -= 1000;
    } while (last == NULL && size >= window);

    if (last == NULL) {
        log_file_printf("check_temp_last: cannot locate the last note.\n");
        return 0;
    }

    wchar_t *a = _wcsdup(note_name);
    wchar_t *b = last;
    int na = (int)wcslen(a);
    int nb = (int)wcslen(b);
    wchar_t *pa = a;
    if (na > 0 && *a != L'.' && *a != L'_') {
        pa = &a[na];
        for (int i = 1; i < na; ++i) {
            if (a[i] == L'.' || a[i] == L'_') {
                pa = &a[i];
                break;
            }
        }
    }
    wchar_t *pb = b;
    if (nb > 0 && *b != L'.' && *b != L'_') {
        pb = &b[nb];
        for (int i = 1; i < nb; ++i) {
            if (b[i] == L'_' || b[i] == L'.') {
                pb = &b[i];
                break;
            }
        }
    }
    *pa = 0;
    *pb = 0;
    log_file_printf("check_temp_last: %ls %ls\n", b, a);
    int ok = wcscmp(b, a) == 0;
    free(a);
    free(b);
    fclose(fp);
    return ok;
}

typedef struct {
    wavtool_args *args;
    int fs;
    wavtool_args **notes;
    float *begins;
    int *cuts;
    int nseg;
    int nsamples;
    float *pcm;
} wavtool_job;

static void wavtool_process_segment(wavtool_job *job)
{
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int chunk = job->nseg / nthreads;
    int rem = job->nseg % nthreads;
    if (tid < rem) {
        ++chunk;
        rem = 0;
    }
    int begin = chunk * tid + rem;
    int end = begin + chunk;
    for (int s = begin; s < end; ++s) {
        int i0 = job->cuts[s];
        int i1 = job->cuts[s + 1] - 1;
        int n = i1 - i0;

        log_tee_printf("main_wavtool: working on segment #%d out of %d\n",
                       s + 1, job->nseg);
        float *rel = (float *)calloc((size_t)(i1 + 1 - i0), sizeof(float));
        float t0 = job->notes[i0]->overlap;
        if (i0 <= i1) {
            float adj = t0 - job->begins[i0];
            for (int k = i0; k <= i1; ++k)
                rel[k - i0] = job->begins[k] + adj;
        }
        float dur = rel[n] + job->notes[i1]->length - job->notes[i1]->overlap;
        log_file_printf(
            "main_wavtool: segment offset = %fsec, segment duration = %fsec\n",
            t0, (double)dur);
        if (i0 > i1) {
            log_tee_printf("main_wavtool: synthesizing...\n");
            log_tee_printf("Error: no note found. Synthesis aborted.\n");
            free(rel);
            continue;
        }

        lmix_data *acc = NULL;
        int nhop = 0;
        llsm_conf saved = {0};
        memset(&saved, 0, sizeof(saved));
        for (int k = i0; k <= i1; ++k) {
            wavtool_args *note = job->notes[k];
            size_t npath = wcslen(note->input_path);
            wchar_t *pack = (wchar_t *)calloc(npath + 10, sizeof(wchar_t));

            wcscpy(pack, note->input_path);
            wcscat(pack, L".llsm.tmp");
            log_file_printf("main_wavtool: reading data pack #%d...\n", k);
            xchg_node *root = xchg_import(pack, 0.0f, 0.0f);
            lmix_data *src = NULL;
            if (root != NULL) {
                src = lmix_from_xchg(xchg_node_get_field(root, "data"));
                if (src != NULL) {
                    xchg_node_destroy(root);
                    log_file_printf("main_wavtool: amplitude enveloping...\n");
                    {
                        lmix_envelope env = wavtool_envelope_from_args(note);
                        lmix_apply_envelope(src, &env);
                    }
                } else {
                    log_tee_printf(
                        "Warning: data pack format not recognized. Skipped. (%ls)\n",
                        pack);
                    xchg_node_destroy(root);
                    free(pack);
                    free_wavtool_args(note);
                    continue;
                }
            } else {
                log_file_printf(
                    "main_wavtool: data pack not found. Searching for raw wav...\n");
                src = lmix_from_file_direct(note->input_path);
                if (src == NULL) {
                    log_tee_printf(
                        "Warning: data pack not found. Skipped. (%ls)\n", pack);
                    free(pack);
                    free_wavtool_args(note);
                    continue;
                }
                src->duration = note->length;
                log_file_printf("main_wavtool: amplitude enveloping...\n");
                {
                    lmix_envelope env = wavtool_envelope_from_args(note);
                    lmix_apply_envelope(src, &env);
                }
                if (src->frames[0]->wave != NULL)
                    src->frames[0]->wave->t0 -= note->offset;
            }

            nhop = src->conf.nhop;
            if (acc == NULL) {
                log_file_printf(
                    "main_wavtool: initializing output data pack...\n");
                int nfrm = (int)(dur * (float)job->fs / (float)nhop) + 4;
                acc = wavtool_empty_lmix(&src->conf, nfrm);
                log_file_printf(
                    "main_wavtool: output data pack initialized.\n");
            } else if (!lmix_conf_compatible(&src->conf, &acc->conf)) {
                log_tee_printf(
                    "Warning: inconsistent model specification. Skipped.\n");
                log_file_printf("main_wavtool: freeing memory...\n");
                lmix_delete(src);
                free(pack);
                free_wavtool_args(note);
                continue;
            }

            log_file_printf(
                "main_wavtool: begin = %f, overlap = %f, duration = %f\n",
                job->begins[k], note->overlap, note->length);
            {
                float inv = 1.0f / (float)nhop;
                int begin = (int)((rel[k - i0] - note->overlap) *
                                  (float)job->fs * inv);
                int overlap = (int)(note->overlap * (float)job->fs * inv);
                int ltruncate = (int)((float)job->fs * note->offset * inv);

                log_tee_printf(
                    "main_wavtool: concatenating data pack #%d...\n", k);
                lmix_superimpose_add(acc, src, begin, overlap, ltruncate, 1);
            }
            log_file_printf("main_wavtool: freeing memory...\n");
            lmix_delete(src);
            free(pack);
            free_wavtool_args(note);
        }

        log_tee_printf("main_wavtool: synthesizing...\n");
        if (acc == NULL) {
            log_tee_printf("Error: no note found. Synthesis aborted.\n");
            free(rel);
            continue;
        } else {
            llsm_parameters param = llsm_init(acc->conf.nnosband);
            param.s_fs = job->fs;
            log_file_printf(
                "main_wavtool: converting LMIX to lv0 LLSM...\n");
            llsm_layer0 *layer0 = lmix_to_layer0(acc);
            log_file_printf(
                "main_wavtool: synthesizing lv0 LLSM to wave...\n");
            llsm_output *out = llsm_layer0_synthesize(param, layer0);
            FP_TYPE *wave = lmix_render_wave(out, acc);
            llsm_deinit(param);
            llsm_delete_layer0(layer0);
            log_file_printf(
                "main_wavtool: appending wave starting from %f sec...\n",
                job->begins[i0] - t0);
            int start = (int)((job->begins[i0] - t0) * (float)job->fs);
            int ny = out->ny;
            for (int i = 0; i < ny; ++i) {
                int pos = start + i;
                if (pos >= 0 && pos < job->nsamples)
                    job->pcm[pos] += wave[i];
            }
            if (log_filename != NULL) {
                log_file_printf(
                    "main_wavtool: outputing LMIX data pack...\n");
                wchar_t suffix[256];
                swprintf(suffix, 256, L".%d", s);
                wchar_t *base = wstr_concat(job->args->output_path, suffix);
                wchar_t *path = (wchar_t *)calloc(wcslen(base) + 10,
                                                  sizeof(wchar_t));
                wcscpy(path, base);
                wcscat(path, L".llsm.tmp");
                xchg_node *node = lmix_to_xchg(acc);
                log_file_printf("output_lmix: xchg_node created.\n");
                xchg_export_counted(node, path);
                log_file_printf("output_lmix: data file exported.\n");
                xchg_node_destroy(node);
                free(path);
                free(base);
            }
            free(wave);
            llsm_delete_output(out);
            lmix_delete(acc);
            free(rel);
        }
    }
}

int main_wavtool(int argc, char **argv)
{
    log_file_printf("main_wavtool: wavtool mode entered.\n");
    wavtool_args *args = parse_wavtool_args(argc, argv);
    log_file_printf("main_wavtool: arguments parsed.\n");
    int fs = g_output_sampling_rate;
    wchar_t *cwd = get_cwd_w();
    wchar_t *in_name = path_get_filename(args->input_path);
    wchar_t *out_name = path_get_filename(args->output_path);

    if (wcscmp(in_name + 2, out_name) != 0) {
        log_file_printf("main_wavtool: chain-calling detected.\n");
        FILE *probe = _wfopen(args->output_path, L"r");
        if (probe != NULL) {
            fclose(probe);
            log_file_printf("main_wavtool: data indexing file detected.\n");
        } else {
            log_file_printf("main_wavtool: entering mode 3...\n");
            probe = _wfopen(args->output_path, L"w");
            fclose(probe);
        }
        log_file_printf("main_wavtool: entering mode 4...\n");
        wavtool_append_index(args->output_path, args);
        log_tee_printf("main_wavtool: data index appended.\n");
        wchar_t *bat = (wchar_t *)calloc(0x400, sizeof(wchar_t));
        wcscpy(bat, cwd);
        wcscat(bat, L"\\temp.bat");
        int last = args->last_note || check_temp_last(bat, in_name);
        if (last) {
            log_file_printf("main_wavtool: entering mode 5...\n");
            if (g_multithread_synthesis == 2)
                omp_set_num_threads(omp_get_max_threads());
            else
                omp_set_num_threads(1);
            int nnotes = 0;
            wavtool_args **notes = wavtool_load_index(args->output_path, &nnotes);
            if (notes == NULL) {
                log_file_printf(
                    "main_wavtool: cannot open indexing file. mode 5 terminated.\n");
                return 1;
            }
            log_file_printf("main_wavtool: data indexing file loaded.\n");
            float *begins = (float *)calloc((size_t)nnotes, sizeof(float));
            int *cuts = (int *)calloc((size_t)nnotes + 2, sizeof(int));
            cuts[0] = 0;
            int nseg = 1;
            float t = 0.0f;
            for (int i = 0; i < nnotes; ++i) {
                wavtool_args *n = notes[i];

                if (i == 0)
                    t = t + n->overlap;
                begins[i] = t;
                t = t + n->length - n->overlap;
                wchar_t *base = path_get_filename(n->input_path);
                log_file_printf("main_wavtool: output & input: %ls, %ls\n",
                                n->output_path, n->input_path);
                log_file_printf("\noffset: %f, length: %f, overlap: %f\n",
                                n->offset, n->length, n->overlap);
                log_file_printf("p1: %f, p2: %f, p3: %f, p4: %f, p5: %f\n",
                                n->p1, n->p2, n->p3, n->p4, n->p5);
                log_file_printf("v1: %f, v2: %f, v3: %f, v4: %f, v5: %f\n",
                                n->v1, n->v2, n->v3, n->v4, n->v5);
                wchar_t *tmp = wstr_concat(n->input_path, L".llsm.tmp");
                if (wcscmp(base, L"R.wav") == 0) {
                    if (i != 0 && nnotes - 1 != i &&
                        n->length > notes[i + 1]->overlap) {
                        log_file_printf("Segmentation occur at #%d\n", i);
                        cuts[nseg++] = i;
                    }
                } else {
                    FILE *tf = _wfopen(tmp, L"r");
                    if (tf == NULL) {
                        if (i != 0 && nnotes - 1 != i &&
                            n->length > notes[i + 1]->overlap) {
                            log_file_printf("Segmentation occur at #%d\n", i);
                            cuts[nseg++] = i;
                        }
                    } else {
                        fclose(tf);
                    }
                }
                free(tmp);
            }
            cuts[nseg] = nnotes;
            int nsamples = (int)((float)fs * t);
            float *pcm = (float *)calloc((size_t)nsamples, sizeof(float));
            log_file_printf(
                "main_wavtool: total duration = %f sec (%d samples)\n",
                t, nsamples);
            wavtool_job job;
            job.args = args;
            job.fs = fs;
            job.notes = notes;
            job.begins = begins;
            job.cuts = cuts;
            job.nseg = nseg;
            job.nsamples = nsamples;
            job.pcm = pcm;
#pragma omp parallel
            {
                wavtool_process_segment(&job);
            }
            log_tee_printf("main_wavtool: outputing wav file...\n");
            char *out_mbs = wchar_to_mbs(args->output_path);
            wavwrite(pcm, nsamples, fs, g_output_bit_depth, out_mbs);
            free(out_mbs);
            free(pcm);
            free(notes);
            free(begins);
            free(cuts);
        }
        free(bat);
    } else {
        log_file_printf("main_wavtool: presamp-style calling detected.\n");
        wchar_t *in_tmp = (wchar_t *)calloc(wcslen(args->input_path) + 10,
                                            sizeof(wchar_t));
        wcscpy(in_tmp, args->input_path);
        wcscat(in_tmp, L".llsm.tmp");
        wchar_t *out_tmp = (wchar_t *)calloc(wcslen(args->output_path) + 10,
                                             sizeof(wchar_t));
        wcscpy(out_tmp, args->output_path);
        wcscat(out_tmp, L".llsm.tmp");
        xchg_node *root = xchg_import(in_tmp, 0.0f, 0.0f);
        if (root == NULL) {
            log_file_printf("Error: cannot load file %ls\n", in_tmp);
            return 1;
        }
        lmix_data *src = lmix_from_xchg(xchg_node_get_field(root, "data"));
        if (src == NULL) {
            log_file_printf(
                "Error: file format not recognized for %ls\n", in_tmp);
            return 1;
        }
        xchg_node_destroy(root);
        {
            lmix_envelope env = wavtool_envelope_from_args(args);
            lmix_apply_envelope(src, &env);
        }
        {
            FILE *exist = _wfopen(args->output_path, L"r");
            if (exist == NULL) {
                int hop = src->conf.nhop;

                log_file_printf("main_wavtool: entering mode 1...\n");
                lmix_data *dst = wavtool_empty_lmix(&src->conf, src->conf.nfrm);
                int ltruncate = (int)((float)fs * args->offset / (float)hop);
                lmix_superimpose_add(dst, src, 0, 0, ltruncate, 1);
                dst->duration = args->length;
                xchg_node *node = lmix_to_xchg(dst);
                xchg_export_counted(node, out_tmp);
                lmix_delete(dst);
                xchg_node_destroy(node);
                FILE *dummy = _wfopen(args->output_path, L"w");
                fwrite("This is a dummy file created just to fool UTAU's cache mechanism "
                       "and you really shouldn't expect to see anything here. The vocals "
                       "are stored in llsm files, if nothing went wrong.\n",
                       1, 0xB3u, dummy);
                fclose(dummy);
            } else {
                fclose(exist);
                log_file_printf("main_wavtool: entering mode 2...\n");
                xchg_node *oroot = xchg_import(out_tmp, 0.0f, 0.0f);
                if (oroot == NULL) {
                    log_file_printf("Error: cannot load file %ls\n", out_tmp);
                    return 1;
                }
                lmix_data *old = lmix_from_xchg(xchg_node_get_field(oroot, "data"));
                if (old == NULL) {
                    log_file_printf(
                        "Error: file format not recognized %ls\n", out_tmp);
                    return 1;
                }
                xchg_node_destroy(oroot);
                {
                    float old_dur = old->duration;
                    float overlap = args->overlap;
                    float new_dur = args->length - overlap + old_dur;
                    int hop = old->conf.nhop;
                    float inv = 1.0f / (float)hop;
                    int nfrm = (int)(inv * ((float)fs * new_dur));
                    int begin = (int)((old_dur - overlap) * (float)fs * inv);
                    int ov = (int)(overlap * (float)fs * inv);
                    int ltruncate = (int)(args->offset * (float)fs * inv);
                    lmix_data *dst = wavtool_empty_lmix(&old->conf, nfrm);

                    lmix_superimpose_add(dst, old, 0, 0, 0, 1);
                    lmix_superimpose_add(dst, src, begin, ov, ltruncate, 1);
                    dst->duration = new_dur;
                    log_tee_printf("main_wavtool: outputing LMIX data pack...\n");
                    xchg_node *node = lmix_to_xchg(dst);
                    xchg_export_counted(node, out_tmp);
                    lmix_delete(old);
                    lmix_delete(dst);
                    xchg_node_destroy(node);
                }
            }
        }
        free(in_tmp);
        free(out_tmp);
        lmix_delete(src);
    }

    log_file_printf("main_wavtool: freeing memory...\n");
    free(cwd);
    free_wavtool_args(args);
    log_file_printf("main_wavtool: wavtool session completed.\n");
    return 0;
}
