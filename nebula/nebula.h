#ifndef MORESAMPLER_NEBULA_H
#define MORESAMPLER_NEBULA_H

// Nebula F0 estimator used by moresampler.
// Algorithm: Kanru Hua, Interspeech 2018 / arXiv:1710.11317.
// nebula_model.c contains embedded GMMs and likelihood-map processing.
// nebula_features.c extracts IF/SNR features.
// nebula_f0.c performs F0 estimation, voicing and refinement.
// Model tables are stored in model/nebula_gmm.inc and model/nebula_lcal.inc.

#include <stddef.h>

typedef struct neb_gmm {
    int nmix;
    int nord;
    float *w;
    float **mu;
    float **sigma;
    float *det11;
    float **inv11;
    float **s21s11i;
    float *condsigma;
    int ready;
} neb_gmm;

// Shared implementation entry points. nebula_est owns the model array and
// likelihood-map allocation for the duration of the estimation.
neb_gmm *neb_gmm_load_band(int band);
void neb_gmm_free(neb_gmm *g);
void neb_make_lmap(
    neb_gmm **bands,
    float **SNR1, float **SNR2, float **IF1, float **IF2,
    int nfrm,
    const float *f_axis,
    float **lmap_rows);

// Result object returned by nebula_est.
// nebula_result_free releases the F0 contour, likelihood rows and result.
typedef struct nebula_result {
    float *f0;              // F0 contour, 0 = unvoiced
    float **lmap_rows;      // per-frame likelihood / feature rows
    void *aux;              // auxiliary data
    float *f_axis;          // candidate F0 axis (Hz or bin centers)
    int nfrm;               // number of frames
    int lmap_width;         // columns in each likelihood row
} nebula_result;

void nebula_result_free(nebula_result *r);

// Estimate F0 and voicing likelihoods from a waveform.
// Processing includes resampling, preprocessing, IF/SNR extraction, GMM
// likelihoods, Viterbi smoothing and F0 refinement.
nebula_result *nebula_est(
    float fmin, float fmax,
    float thop,
    int use_peak_likelihood,
    int skip_voicing,
    float dither_level,
    float voicing_time_const,
    float viterbi_sigma,
    float lmap_smooth,
    float voicing_mu,
    float voicing_var,
    int unused12, int unused13, int unused14,
    float *wave, int nx,
    float fs);

// Estimate an F0 contour from a waveform using the LLSM analysis parameters.
// Returns a newly allocated array of length *out_nfrm; the caller frees it.
float *f0_from_wave(
    float fmin, float fmax,
    float thop,
    int use_peak_likelihood,
    int skip_voicing,
    float dither_level,
    float voicing_time_const,
    float viterbi_sigma,
    float lmap_smooth,
    float voicing_mu,
    float voicing_var,
    int n_pass, int nhar, int kalman_win,
    float *x, int nx,
    int fs,
    int *out_nfrm);

// Apply optional post-refinement to an F0 contour.
// Leading analysis parameters match f0_from_wave / nebula_est; only
// thop, n_pass, nhar, kalman_win and the waveform/F0 inputs are used.
float *biased_f0_estimation(
    float fmin, float fmax,
    float thop,
    int use_peak_likelihood,
    int skip_voicing,
    float dither_level,
    float voicing_time_const,
    float viterbi_sigma,
    float lmap_smooth,
    float voicing_mu,
    float voicing_var,
    int n_pass, int nhar, int kalman_win,
    float *wave, int nx,
    float fs,
    float *f0_in, int nfrm);

// Extract IF/SNR features for each analysis band.
void estimate_ifsnr(
    float *x, int nx,
    const int *naxis, int nfrm,
    const float *fc, const float *fres, int nband,
    const char *wtype, const char *stype,
    float **SNR1, float **IF1,
    float **SNR2, float **IF2);

// Compute voicing posterior probabilities.
float *nebula_voicing_posterior(
    int Count, int use_uv_offset, const float *x,
    float mu, float var, float *out_logZ);

#endif /* MORESAMPLER_NEBULA_H */
