#ifndef MORESAMPLER_RESAMPLER_H
#define MORESAMPLER_RESAMPLER_H

#include <wchar.h>

// resampler_args holds everything parsed from UTAU's resampler command line:
// moresampler <in.wav> <out.wav> <pitch> <velocity> [flags] [offset]
// [length] [consonant] [cutoff] [volume] [modulation] [!tempo] [pitchbend]
// plus the fully-decoded, range-clamped flag values (g, t, Mt, ...).
typedef struct resampler_args {
    wchar_t *in_path;           // input wav path (wchar_t buffer)
    wchar_t *out_path;          // output path (wchar_t buffer)
    char *pitch;                // note name / pitch string (strdup'd)
    float velocity;             // consonant velocity

    char *flags;                // raw flag sequence string

    float offset;               // region offset (sec)
    float length_require;       // required length (sec), -1 if unset
    float fixed_length;         // consonant/fixed length (sec)
    float end_blank;            // cutoff / end blank (sec; <0 = length)
    float volume;               // volume (0..100+, default 100)
    float modulation;           // modulation (default 0)
    float tempo;                // "!"-prefixed BPM for pitchbend timing
    float unknow12;

    float *pitch_bend;          // decoded pitch-bend curve (cents)
    int pitch_bend_size;

    // Decoded & clamped flag values. Indices match flag_parse name table.
    float g;                    // gender          [-100, 100]
    float t;                    // pitch offset    [-1200, 1200]
    float Y;                    // (raw)
    int e;                      // +1 if "e", -1 if "Me", else 0
    float b;                    // unvoiced gain   [~-20, 100]
    float A;                    // amplitude mod   [-100, 100]
    float P;                    // compressor      [0, 100]
    int u;                      // "u" flag set
    float Mt;                   // tension         [-100, 100]
    float Mb;                   // breathiness     [-100, 100]
    float Md;                   // dryness         [-100, 100]
    float Mm;                   // model interp    [0, 100]
    int Ms;                     // stabilization   [0, 10] after rint
    float ME;                   // formant emphasis[-100, 100]
    float Mo;                   // openness        [-100, 100]
    float Mr;                   // resonance       [-100, 100]
    float MC;                   // coarseness      [0, 100]
    float MG;                   // growl           [0, 100]
    float MD;                   // distortion      [0, 100]
    float MP;                   // (Mp) peak       [0, 100]
    float _s;                   // (raw)

    int error_code;             // 1 => analysis-only (output == L"nul")
    int skip_llsm;              // skip .llsm cache when "_f" is present
} resampler_args;

resampler_args *parse_resampler_args(int argc, char **argv);
void free_resampler_args(resampler_args *args);

int main_resampler(int argc, char **argv);

#endif // MORESAMPLER_RESAMPLER_H
