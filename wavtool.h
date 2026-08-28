#ifndef MORESAMPLER_WAVTOOL_H
#define MORESAMPLER_WAVTOOL_H

#include <wchar.h>

// wavtool_args stores paths, timing, envelope points, and the last-note flag.
typedef struct {
    wchar_t *output_path;
    wchar_t *input_path;
    float offset;
    float length;
    float overlap;
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
    int last_note;
} wavtool_args;

// Parse wavtool command-line arguments.
wavtool_args *parse_wavtool_args(int argc, char **argv);
// Release wavtool arguments and their paths.
void free_wavtool_args(wavtool_args *args);
// Run wavtool processing.
int main_wavtool(int argc, char **argv);

#endif // MORESAMPLER_WAVTOOL_H
