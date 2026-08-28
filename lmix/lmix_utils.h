#ifndef MORESAMPLER_LMIX_UTILS_H
#define MORESAMPLER_LMIX_UTILS_H

#define LIBAIFF_NOCOMPAT 1


#include <wchar.h>

int sign_change_span(FP_TYPE *x, int n, int i);
int hop_from_fs(int fs);

int midi_from_note(const char *note);
FP_TYPE midi_to_hz(int midi);
FP_TYPE hz_to_mel(FP_TYPE hz);
int parse_note_f0(const char *note, FP_TYPE *f0_out);
FP_TYPE f0_from_dirname(const wchar_t *path);

size_t path_dir_prefix_len(const wchar_t *path);
time_t file_mtime(const wchar_t *path);
int frq_index(int sample, FP_TYPE frq_hop, int count);

FP_TYPE lerp_phase(FP_TYPE a, FP_TYPE b, FP_TYPE t);

// Read .wav via wavread or AIFF via libaiff; multi-channel AIFF keeps channel 0.
FP_TYPE *audioread(
    const wchar_t *path, int *fs, int *nbit, int *nsample);

#endif /* MORESAMPLER_LMIX_UTILS_H */
