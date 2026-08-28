#ifndef MORESAMPLER_CONFIG_H
#define MORESAMPLER_CONFIG_H

#include <wchar.h>

// Runtime configuration globals; defaults match the data section.
extern int g_multithread_synthesis; // 0=off, 1=on, 2=full

extern int g_output_sampling_rate;
extern int g_output_bit_depth;
extern int g_resampler_compatibility;
extern int g_auto_update_llsm_mrq;
extern int g_synthesis_utau_style_normalization;
extern int g_synthesis_loudness_preservation;
extern int g_synthesis_duration_extension_method;
extern int g_load_frq;
extern int g_analysis_noise_reduction;
extern int g_analysis_anti_distortion;
extern int g_analysis_suppress_subharmonics;
extern int g_analysis_biased_f0_estimation;
extern int g_analysis_f0_range_from_path;
extern float g_analysis_f0_min;
extern float g_analysis_f0_max;
extern float g_analysis_llsm_max_voiced_frequency;
extern float g_analysis_llsm_max_noise_frequency;
extern float g_analysis_llsm_mbn_band0;
extern float g_analysis_llsm_mbn_band1;
extern float g_analysis_llsm_mbn_band2;
extern int g_analysis_llsm_max_num_harmonics;
extern int g_analysis_llsm_energy_envelope_num_harmonics;

// is_global=1 loads moreconfig.txt at startup; 0 loads a per-voice override.
int load_config(const wchar_t *path, int is_global);

#endif // MORESAMPLER_CONFIG_H
