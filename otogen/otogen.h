#ifndef MORESAMPLER_OTOGEN_H
#define MORESAMPLER_OTOGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wchar.h>
#include <windows.h>

#include "../external/ciglet/ciglet.h"
#include "../external/liblrhsmm/common.h"
#include "../external/liblrhsmm/data.h"
#include "../external/liblrhsmm/model.h"
#include "../external/liblrhsmm/estimate.h"
#include "../external/liblrhsmm/inference.h"
#include "../external/cJSON/cJSON.h"

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

typedef struct otogen_hnode {
    wchar_t *key;
    int index;
    struct otogen_hnode *next;
} otogen_hnode;

typedef struct {
    int count;
    int nbucket;
    otogen_hnode **head;
    otogen_hnode **tail;
} otogen_hash;

typedef struct otogen_plist {
    wchar_t *phoneme;
    int unused;
    struct otogen_plist *next;
} otogen_plist;

typedef struct {
    otogen_hash *map;
    otogen_plist **chains;
    size_t max_glen;
} otogen_g2p;

typedef struct {
    otogen_hash *all;
    otogen_hash *consonant;
    otogen_hash *vowel;
    otogen_hash *voiced;
    otogen_hash *unvoiced;
    otogen_hash *semi;
    otogen_hash *silence;
    wchar_t **names;
} otogen_phdict;

typedef struct {
    FP_TYPE *tbegin;
    FP_TYPE *tend;
    wchar_t **phoneme;
    wchar_t **grapheme;
    int nseg;
} otogen_label;

typedef struct {
    float *mean;
    int ndim;
    float weight;
} otogen_bmodel_cell;

typedef struct {
    otogen_bmodel_cell ***cell;
    int nclass;
    int ndim;
} otogen_bmodel;

typedef struct {
    int nstate;
    int nstream;
    int ncep;
    int nmel;
    int nmel_delta;
    FP_TYPE fs;
    FP_TYPE fmin;
    FP_TYPE fmax;
    FP_TYPE hop_sec;
    FP_TYPE win_sec;
    FP_TYPE fft_sec;
    FP_TYPE silence_thresh;
    FP_TYPE pad_sec;
    FP_TYPE silence_floor;
    int nskip;
    otogen_g2p *g2p;
    otogen_phdict *phdict;
    int nfile;
    wchar_t **wav_paths;
    wchar_t **aliases;
    otogen_label **labels;
    int *infer_nt;
    FP_TYPE *trim;
    lrh_segset *segset;
    lrh_observset *observset;
    FP_TYPE ***band_e;
    lrh_model *model;
    otogen_bmodel *boundary_model;
} otogen_ctx;

int main_oto_generation(int argc, char **argv);

unsigned otogen_hash_wchar(const wchar_t *s, int nbucket);
otogen_hash *otogen_hash_create(int nbucket);
void otogen_hash_insert(otogen_hash *h, const wchar_t *key, int index);
otogen_hnode *otogen_hash_find(otogen_hash *h, const wchar_t *key);
int otogen_hash_has(otogen_hash *h, const wchar_t *key);
void otogen_hash_free(otogen_hash *h);
otogen_hash *otogen_hash_from_cjson_array(cJSON *arr, otogen_hash *all);

otogen_g2p *otogen_g2p_from_csv(const wchar_t *text);
void otogen_g2p_free(otogen_g2p *g);
otogen_label *otogen_g2p_apply(otogen_g2p *g, const wchar_t *text);
void otogen_label_free(otogen_label *lab);

otogen_phdict *otogen_phdict_from_json(const char *json);
void otogen_phdict_free(otogen_phdict *d);
int otogen_is_consonant(otogen_phdict *d, const wchar_t *ph);
int otogen_is_vowel(otogen_phdict *d, const wchar_t *ph);
int otogen_is_voiced(otogen_phdict *d, const wchar_t *ph);
int otogen_is_unvoiced(otogen_phdict *d, const wchar_t *ph);
int otogen_is_semi(otogen_phdict *d, const wchar_t *ph);
int otogen_is_silence(otogen_phdict *d, const wchar_t *ph);

otogen_ctx *otogen_ctx_create(void);
void otogen_ctx_free(otogen_ctx *ctx);
int otogen_scan_wavs(otogen_ctx *ctx, const wchar_t *dir);
int otogen_aliases_from_wavs(otogen_ctx *ctx);
int otogen_analyze_audio(otogen_ctx *ctx);
void otogen_init_segmentation(otogen_ctx *ctx);
void otogen_init_model(otogen_ctx *ctx, int skip_stats);
otogen_bmodel *otogen_bmodel_from_table(void);
void otogen_bmodel_free(otogen_bmodel *m);
void otogen_pool_vowel_gmms(otogen_ctx *ctx);
void otogen_train(otogen_ctx *ctx, int niter, int hsmm);
void otogen_edit_silence_floor(otogen_ctx *ctx);
void otogen_align(otogen_ctx *ctx, int geometric);
void otogen_refine_alignment(otogen_ctx *ctx, int mode);
void otogen_fill_label_times(otogen_ctx *ctx);

wchar_t *otogen_fastlz_to_wchar(const void *src, size_t clen, int dlen);
const char *otogen_luac_reader(lua_State *L, void *ud, size_t *sz);

int otogen_lua_printstr(lua_State *L);
int otogen_lua_readstr(lua_State *L);
int otogen_lua_is_consonant(lua_State *L);

#endif // MORESAMPLER_OTOGEN_H
