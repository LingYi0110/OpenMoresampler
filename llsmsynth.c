#include "llsmsynth.h"

#include <stdlib.h>

#include "config.h"
#include "external/ciglet/ciglet.h"
#include "external/libllsm/llsm.h"
#include "utils.h"
#include "lmix.h"
#include "xchg.h"

// Load an LMIX pack, synthesize it, and write output.wav.
int main_llsmsynth(int argc, char **argv)
{
    int fs = g_output_sampling_rate;
    wchar_t *path = get_argv_w(argc, argv, 1);
    xchg_node *root = xchg_import(path, 0.0f, 0.0f);
    lmix_data *lmix = lmix_from_xchg(xchg_node_get_field(root, "data"));
    xchg_node_destroy(root);

    llsm_parameters param = llsm_init(lmix->conf.nnosband);
    llsm_layer0 *layer0 = lmix_to_layer0(lmix);
    param.s_fs = fs;
    llsm_output *out = llsm_layer0_synthesize(param, layer0);
    llsm_deinit(param);

    wavwrite(out->y, out->ny, fs, g_output_bit_depth, "output.wav");

    lmix_delete(lmix);
    llsm_delete_output(out);
    llsm_delete_layer0(layer0);
    free(path);
    return 0;
}
