#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wchar.h>
#include <math.h>

#include "logging.h"
#include "utils.h"
#include "flag.h"


int g_multithread_synthesis = 0;

int g_output_sampling_rate = 44100;
int g_output_bit_depth = 16;
int g_resampler_compatibility = 1;
int g_auto_update_llsm_mrq = 0;
int g_synthesis_utau_style_normalization = 2; // 0=off,1=voiced,2=full
int g_synthesis_loudness_preservation = 1;
int g_synthesis_duration_extension_method = 0; // 0=auto,1=stretch,2=?
int g_load_frq = 0; // 0=off,1=on,2=strict
int g_analysis_noise_reduction = 0;
int g_analysis_anti_distortion = 0;
int g_analysis_suppress_subharmonics = 0;
int g_analysis_biased_f0_estimation = 0;
int g_analysis_f0_range_from_path = 1;
float g_analysis_f0_min = 80.0f;
float g_analysis_f0_max = 800.0f;
float g_analysis_llsm_max_voiced_frequency = 12000.0f;
float g_analysis_llsm_max_noise_frequency = 24000.0f;
float g_analysis_llsm_mbn_band0 = 3000.0f;
float g_analysis_llsm_mbn_band1 = 6000.0f;
float g_analysis_llsm_mbn_band2 = 12000.0f;
int g_analysis_llsm_max_num_harmonics = 80;
int g_analysis_llsm_energy_envelope_num_harmonics = 4;

// Per-key sticky override flags (set when line starts with '*').
static int ov_utau_norm;
static int ov_loudness;
static int ov_duration_ext;
static int ov_load_frq;
static int ov_noise_red;
static int ov_anti_dist;
static int ov_suppress_sub;
static int ov_biased_f0;
static int ov_f0_range;
static int ov_f0_min;
static int ov_f0_max;
static int ov_max_voiced;
static int ov_max_noise;
static int ov_mbn0;
static int ov_mbn1;
static int ov_mbn2;
static int ov_max_harm;
static int ov_env_harm;

static int wcs_on(const wchar_t *v)
{
    return wcsncmp(v, L"on", 2) == 0;
}

static float parse_float_w(const wchar_t *s)
{
    // Use the wide floating-point conversion path.
    return (float)wcstod(s, NULL);
}

// Load global or voice-local configuration values.
int load_config(const wchar_t *path, int is_global)
{
    FILE *fp = _wfopen(path, L"r");
    if (!fp) {
        if (is_global)
            log_tee_printf("load_config: configuration file not found.\n");
        return is_global;
    }

    for (;;) {
        WCHAR *line = fgets_utf8_to_wchar(fp);
        if (!line)
            break;

        size_t len = wcslen(line);
        if (len <= 2 || line[0] == L'#') {
            free(line);
            continue;
        }

        // Split on the first space into key and value.
        WCHAR *p = line;
        while (*p && *p != L' ')
            p++;
        if (*p == L' ')
            *p++ = 0;
        while (*p == L' ')
            p++;

        int sticky = 0;
        WCHAR *key = line;
        if (key[0] == L'*') {
            sticky = 1;
            key++;
        }
        WCHAR *val = p;
        // Trim trailing CR/LF from the value.
        size_t vl = wcslen(val);
        while (vl > 0 && (val[vl - 1] == L'\n' || val[vl - 1] == L'\r')) {
            val[--vl] = 0;
        }

        if (is_global) {
            if (wcscmp(key, L"output-sampling-rate") == 0) {
                g_output_sampling_rate = (int)wcstol(val, NULL, 0);
            } else if (wcscmp(key, L"output-bit-depth") == 0) {
                g_output_bit_depth = (int)wcstol(val, NULL, 0);
            } else if (wcscmp(key, L"resampler-compatibility") == 0) {
                g_resampler_compatibility = wcs_on(val);
            } else if (wcscmp(key, L"dump-log-file") == 0) {
                if (wcsncmp(val, L"off", 3) == 0) {
                    log_filename = NULL;
                } else {
                    log_filename = _wcsdup(val);
                    if (log_filename) {
                        // Strip the trailing newline from the log filename.
                        size_t n = wcslen(log_filename);
                        if (n > 0)
                            log_filename[n - 1] = 0;
                        log_fp = _wfopen(log_filename, L"a");
                        if (!log_fp) {
                            log_fp = _wfopen(log_filename, L"w+");
                            if (log_fp)
                                fseek(log_fp, 0, SEEK_END);
                            else
                                log_tee_printf("Error: cannot create/append log file!\n");
                        }
                    }
                }
            } else if (wcscmp(key, L"multithread-synthesis") == 0) {
                if (wcs_on(val))
                    g_multithread_synthesis = 1;
                else if (wcsncmp(val, L"full", 4) == 0)
                    g_multithread_synthesis = 2;
                else
                    g_multithread_synthesis = 0;
            } else if (wcscmp(key, L"auto-update-llsm-mrq") == 0) {
                g_auto_update_llsm_mrq = wcs_on(val);
            } else if (wcsncmp(key, L"meta-flag-", 10) == 0) {
                int idx = (int)wcstol(key + 10, NULL, 0);
                if ((unsigned)idx > 0xFFF) {
                    log_tee_printf("Error: invalid meta flag number '%ls'.\n", key);
                } else {
                    if (idx >= g_flag_capacity) {
                        int newcap = g_flag_capacity * 2;
                        if (newcap < idx + 1)
                            newcap = idx + 1;
                        g_flag_tables = realloc(g_flag_tables, sizeof(flag_table *) * (size_t)newcap);
                        for (int i = g_flag_capacity; i < newcap; i++)
                            g_flag_tables[i] = NULL;
                        g_flag_capacity = newcap;
                    }
                    if (g_flag_tables[idx] == NULL)
                        g_flag_tables[idx] = flag_create();
                    char *mb = wchar_to_mbs(val);
                    flag_parse(100000.0f, g_flag_tables[idx], mb);
                    free(mb);
                }
            }
        }

        // Apply keys accepted in both modes with sticky override logic.
        if (wcscmp(key, L"synthesis-utau-style-normalization") == 0) {
            if (is_global || !ov_utau_norm) {
                int v;
                if (wcsncmp(val, L"full", 4) == 0)
                    v = 2;
                else
                    v = wcsncmp(val, L"voiced", 6) == 0;
                g_synthesis_utau_style_normalization = v;
                if (sticky)
                    ov_utau_norm = 1;
            }
        } else if (wcscmp(key, L"synthesis-loudness-preservation") == 0) {
            if (is_global || !ov_loudness) {
                g_synthesis_loudness_preservation = wcs_on(val);
                if (sticky)
                    ov_loudness = 1;
            }
        } else if (wcscmp(key, L"synthesis-duration-extension-method") == 0) {
            if (is_global || !ov_duration_ext) {
                if (wcsncmp(val, L"auto", 4) == 0)
                    g_synthesis_duration_extension_method = 0;
                else
                    g_synthesis_duration_extension_method =
                        (wcsncmp(val, L"stretch", 7) != 0) + 1;
                if (sticky)
                    ov_duration_ext = 1;
            }
        } else if (wcscmp(key, L"load-frq") == 0) {
            if (is_global || !ov_load_frq) {
                if (wcs_on(val))
                    g_load_frq = 1;
                else if (wcsncmp(val, L"strict", 6) == 0)
                    g_load_frq = 2;
                else
                    g_load_frq = 0;
                if (sticky)
                    ov_load_frq = 1;
            }
        } else if (wcscmp(key, L"analysis-noise-reduction") == 0) {
            if (is_global || !ov_noise_red) {
                g_analysis_noise_reduction = wcs_on(val);
                if (sticky)
                    ov_noise_red = 1;
            }
        } else if (wcscmp(key, L"analysis-anti-distortion") == 0) {
            if (is_global || !ov_anti_dist) {
                g_analysis_anti_distortion = wcs_on(val);
                if (sticky)
                    ov_anti_dist = 1;
            }
        } else if (wcscmp(key, L"analysis-suppress-subharmonics") == 0) {
            if (is_global || !ov_suppress_sub) {
                g_analysis_suppress_subharmonics = wcs_on(val);
                if (sticky)
                    ov_suppress_sub = 1;
            }
        } else if (wcscmp(key, L"analysis-biased-f0-estimation") == 0) {
            if (is_global || !ov_biased_f0) {
                g_analysis_biased_f0_estimation = wcs_on(val);
                if (sticky)
                    ov_biased_f0 = 1;
            }
        } else if (wcscmp(key, L"analysis-f0-range-from-path") == 0) {
            if (is_global || !ov_f0_range) {
                g_analysis_f0_range_from_path = wcs_on(val);
                if (sticky)
                    ov_f0_range = 1;
            }
        } else if (wcscmp(key, L"analysis-f0-min") == 0) {
            if (is_global || !ov_f0_min) {
                g_analysis_f0_min = parse_float_w(val);
                if (sticky)
                    ov_f0_min = 1;
            }
        } else if (wcscmp(key, L"analysis-f0-max") == 0) {
            if (is_global || !ov_f0_max) {
                g_analysis_f0_max = parse_float_w(val);
                if (sticky)
                    ov_f0_max = 1;
            }
        } else if (wcscmp(key, L"analysis-llsm-max-voiced-frequency") == 0) {
            if (is_global || !ov_max_voiced) {
                g_analysis_llsm_max_voiced_frequency = parse_float_w(val);
                if (sticky)
                    ov_max_voiced = 1;
            }
        } else if (wcscmp(key, L"analysis-llsm-max-noise-frequency") == 0) {
            if (is_global || !ov_max_noise) {
                g_analysis_llsm_max_noise_frequency = parse_float_w(val);
                if (sticky)
                    ov_max_noise = 1;
            }
        } else if (wcscmp(key, L"analysis-llsm-mbn-band0") == 0) {
            if (is_global || !ov_mbn0) {
                g_analysis_llsm_mbn_band0 = parse_float_w(val);
                if (sticky)
                    ov_mbn0 = 1;
            }
        } else if (wcscmp(key, L"analysis-llsm-mbn-band1") == 0) {
            if (is_global || !ov_mbn1) {
                g_analysis_llsm_mbn_band1 = parse_float_w(val);
                if (sticky)
                    ov_mbn1 = 1;
            }
        } else if (wcscmp(key, L"analysis-llsm-mbn-band2") == 0) {
            if (is_global || !ov_mbn2) {
                g_analysis_llsm_mbn_band2 = parse_float_w(val);
                if (sticky)
                    ov_mbn2 = 1;
            }
        } else if (wcscmp(key, L"analysis-llsm-max-num-harmonics") == 0) {
            if (is_global || !ov_max_harm) {
                g_analysis_llsm_max_num_harmonics = (int)wcstol(val, NULL, 0);
                if (sticky)
                    ov_max_harm = 1;
            }
        } else if (wcscmp(key, L"analysis-llsm-energy-envelope-num-harmonics") == 0) {
            if (is_global || !ov_env_harm) {
                g_analysis_llsm_energy_envelope_num_harmonics =
                    (int)wcstol(val, NULL, 0);
                if (sticky)
                    ov_env_harm = 1;
            }
        }

        free(line);
    }

    fclose(fp);
    return 0;
}
