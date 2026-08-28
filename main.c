#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <wchar.h>

#include "logging.h"
#include "utils.h"
#include "config.h"
#include "flag.h"
#include "resampler.h"
#include "llsmsynth.h"
#include "wavtool.h"
#include "otogen/otogen.h"


#define MODE_LLSMSYNTH 0
#define MODE_RESAMPLER 10
#define MODE_WAVTOOL 11
#define MODE_OTO_GENERATION 12
#define MODE_USAGE 13


static int detect_mode(int argc, char **argv)
{
    if (argc <= 1)
        return MODE_USAGE;

    if (argc == 2) {
        wchar_t *wpath = mbs_to_wchar(argv[1]);
        int is_dir = path_is_directory(wpath);
        free(wpath);
        return is_dir ? MODE_OTO_GENERATION : MODE_LLSMSYNTH;
    }

    const char *arg4 = argv[4];
    size_t len = strlen(arg4);
    int marker = 0;
    for (size_t i = 0; i < len; i++) {
        if (arg4[i] == '@') {
            marker = 1;
            break;
        }
    }
    if ((unsigned char)(argv[3][0] - '0') <= 9)
        marker = 1;

    if (strcmp(argv[2], "nul") == 0 || marker == 0)
        return MODE_RESAMPLER;
    return MODE_WAVTOOL;
}

int main(int argc, char **argv)
{
    flag_init();

    printf("Moresampler %d.%d.%d - More than a resampler.\n", 0, 8, 4);
    puts("Copyright Kanru Hua, 2015-2018. All rights reserved.");

    // Load moreconfig.txt next to the executable
    wchar_t *exe_path = get_argv_w(argc, argv, 0);
    *path_get_filename(exe_path) = L'\0';
    wchar_t *cfg_path = wstr_concat(exe_path, L"moreconfig.txt");
    load_config(cfg_path, 1);
    free(exe_path);
    free(cfg_path);

    time_t t = time(NULL);
    char *time_char = asctime(gmtime(&t));
    log_file_printf("main: new session started on %s", time_char);

     // If a log file path is set, force single-thread synthesis.
     // Otherwise, if multithread is off, also pin to 1 thread.
    if (log_filename != NULL) {
        g_multithread_synthesis = 0;
        omp_set_num_threads(1);
    } else if (g_multithread_synthesis == 0) {
        omp_set_num_threads(1);
    }

    int mode = detect_mode(argc, argv);
    int rc;

    if (mode == MODE_WAVTOOL) {
        rc = main_wavtool(argc, argv);
    } else if (mode == MODE_LLSMSYNTH) {
        rc = main_llsmsynth(argc, argv);
    } else if (mode == MODE_RESAMPLER) {
        rc = main_resampler(argc, argv);
    } else if (mode == MODE_OTO_GENERATION) {
        rc = main_oto_generation(argc, argv);
    } else {
        rc = 0;
        log_tee_printf(
            "Moresampler only works as a backend for UTAU. "
            "Please run Moresampler from the host software.");
        getchar();
    }

    log_close();
    return rc;
}
