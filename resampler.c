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
#include "logging.h"
#include "utils.h"
#include "lmix.h"
#include "xchg.h"

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
        lmix = lmix_from_file(
            args->in_path, fmax(args->offset - 0.05, 0.0),
            args->end_blank < 0.0f ? args->offset - args->end_blank : -args->end_blank,
            args->skip_llsm);
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

    // Automatic mapping does not consume the optional explicit-map pointers.
    lmix_time_mapping mapping;
    mapping.nmap = 0;
    mapping.map_src = NULL;
    mapping.map_dst = NULL;
    float total_dur = (float)lmix->conf.nfrm * lmix->conf.thop;
    float offset = args->offset;
    float stretch_end;
    if (args->end_blank < 0.0f)
        stretch_end = args->end_blank - offset + total_dur;
    else
        stretch_end = args->end_blank;
    float fixed_length = args->fixed_length;
    mapping.stretch = pow(2.0, 1.0 - 0.01 * args->velocity) * fixed_length;
    float length = args->length_require;
    if (length < 0.0f)
        length = total_dur;

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
            int skip = (int)round((double)wave->fs * (offset - 0.01));
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
        mapping.offset = offset;
        mapping.stretch_end = stretch_end;
        mapping.fixed_length = fixed_length;
        mapping.total_dur = total_dur;
        mapping.length = length;
        lmix2 = lmix_modify_data(lmix, args, mapping);
    }

    if (g_synthesis_utau_style_normalization != 0) {
        log_file_printf("main_resampler: normalizing...\n");
        lmix_normalize_data(lmix2, 0.005 * args->volume, args->P);
    } else {
        lmix_scale_volume(lmix2, 0.01 * args->volume);
    }

    if (args->Mb > 0.0f || args->MC > 0.0f) {
        log_file_printf("main_resampler: applying 'Mb' and/or 'MC' flag...\n");
        float mbmc = args->MC;
        if (mbmc < args->Mb)
            mbmc = args->Mb;
        double gain = log(1.0000000001 - mbmc * 0.01);
        for (int i = 0; i < lmix2->conf.nfrm; i++) {
            lmix_frame *frame = lmix2->frames[i];
            int n = frame->nhar * frame->vt_stride;
            for (int k = 0; k < n; k++)
                frame->vt_magn[k] += gain;
        }
    }

    if (g_resampler_compatibility != 0) {
        llsm_parameters param = llsm_init(lmix2->conf.nnosband);
        int fs = g_output_sampling_rate;
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
        const wchar_t *out_path = args->out_path;
        wchar_t *tmp_path = calloc(wcslen(out_path) + 10, sizeof(wchar_t));
        wcscpy(tmp_path, out_path);
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
