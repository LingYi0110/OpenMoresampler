#ifndef MORESAMPLER_LMIX_LF_H
#define MORESAMPLER_LMIX_LF_H

#include "../external/ciglet/ciglet.h"

// LF spectrum packed params {T0,Te,Tp,Ta,wg,sin,cos,e,A,a,E0}
typedef struct {
    FP_TYPE T0;
    FP_TYPE Te;
    FP_TYPE Tp;
    FP_TYPE Ta;
    FP_TYPE wg;
    FP_TYPE sin_wgTe;
    FP_TYPE cos_wgTe;
    FP_TYPE e;
    FP_TYPE A;
    FP_TYPE a;
    FP_TYPE E0;
} lmix_lfparam;

lmix_lfparam lmix_lfparam_from_model(lfmodel model);
FP_TYPE *lmix_lfmodel_spectrum(
    lfmodel model, FP_TYPE *freq, int nf, FP_TYPE *dst_phase);

#endif
