#ifndef MORESAMPLER_LMIX_H
#define MORESAMPLER_LMIX_H


#include "llsm.h"
#include "xchg.h"


typedef struct lmix_frame {
    FP_TYPE f0;                // fundamental frequency (Hz)
    int nhar;                    // number of harmonics
    int vt_stride;              // vocal-tract magnitude samples stored per harmonic
    FP_TYPE rd;                 // Rd parameter of the LF glottal model
    FP_TYPE *vs_har_phse;         // voice-source harmonic phases
    FP_TYPE *vt_magn;        // vocal-tract log-magnitude, [nhar * vt_stride]
    FP_TYPE vt_freq_step;       // f0 / vt_stride (Hz)
    llsm_nosframe *noise;       // noise/stochastic component
    llsm_wave *wave;            // optional time-domain PCM
    FP_TYPE distortion;         // MD / 100
    FP_TYPE growl;              // MG / 100
} lmix_frame;

typedef struct lmix_data {
    int version;
    FP_TYPE duration;
    llsm_conf conf;
    lmix_frame **frames;
} lmix_data;

extern int g_lmix_nband_default;
extern int g_lmix_nfft;

lmix_data *lmix_create(llsm_conf conf);

lmix_frame *lmix_create_frame(llsm_conf conf, int nfft, int nhar);

void lmix_copy_frame(lmix_frame *dst, lmix_frame *src);
void lmix_frame_delete(lmix_frame *frame);
void lmix_delete(lmix_data *lmix);

lmix_data *lmix_from_file(
    const wchar_t *path, float start, float end,
    int skip_llsm);
lmix_data *lmix_from_file_direct(
    const wchar_t *path);
lmix_data *lmix_from_wave(
    FP_TYPE *samples, int nsample, int fs,
    FP_TYPE *f0, int nf0, int nhop);

llsm_layer0 *lmix_to_layer0(lmix_data *lmix);
xchg_node *lmix_to_xchg(lmix_data *lmix);
lmix_data *lmix_from_xchg(xchg_node *node);

// nmap/map_src/map_dst are mapping parameters; moresampler currently passes 0, NULL, NULL.
struct resampler_args;
lmix_data *lmix_modify_data(
    lmix_data *lmix, struct resampler_args *args,
    float offset, float stretch_end, float fixed_length, float total_dur,
    float stretch, float length,
    int nmap, const void *map_src, const void *map_dst);

void lmix_unwrap_phase_track(
    FP_TYPE *phase, FP_TYPE *weight, int n, int nbin, int ntran);

/*
 * UTAU wavtool envelope.  Through v5 its layout matches wavtool_args:
 * two path pointers, offset, length, overlap, p1..p5, then v1..v5.
 * The point order is p1 -> p2 -> p5 ... p3 <- p4, not numeric order.
 */
typedef struct lmix_envelope {
    void *unused_output_path;
    void *unused_input_path;
    float t0;
    float duration;
    float unused_overlap;
    float p1;
    float p2;
    float p3;
    float p4;
    float p5;
    float v1;
    float v2;
    float v3;
    float v4;
    float v5;
} lmix_envelope;

void lmix_apply_envelope(lmix_data *lmix, const lmix_envelope *env);
void lmix_scale_volume(lmix_data *lmix, float scale);
void lmix_normalize_data(lmix_data *lmix, float scale, float P);
void lmix_superimpose_add(
    lmix_data *dst, lmix_data *src,
    int begin, int overlap, int ltruncate, int flag);
FP_TYPE *lmix_render_wave(llsm_output *out, lmix_data *lmix);

#endif /* MORESAMPLER_LMIX_H */
