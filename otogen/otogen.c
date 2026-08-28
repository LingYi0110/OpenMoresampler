#include "otogen.h"
#include "otogen_resources.inc"

#include "../external/FastLZ/fastlz.h"
#include "../logging.h"
#include "../utils.h"

static int g_otogen_luac_fed;
static wchar_t *g_g2p_hiragana;
static wchar_t *g_g2p_romaji;
static wchar_t *g_g2p_english;


static const luaL_Reg otogen_ext_regs[] = {
    {"is_consonant", otogen_lua_is_consonant},
    {"printstr", otogen_lua_printstr},
    {"readstr", otogen_lua_readstr},
    {NULL, NULL}
};

int otogen_lua_printstr(lua_State *L)
{
    const char *s = lua_tolstring(L, 1, NULL);
    wchar_t *w = utf8_to_wchar(s);
    char *mbs = wchar_to_mbs(w);
    log_tee_printf("%s", mbs);
    free(mbs);
    free(w);
    lua_settop(L, -2);
    return 1;
}


int otogen_lua_readstr(lua_State *L)
{
    char buf[256];
    memset(buf, 0, 256);
    fgets(buf, 256, stdin);
    int i = 255;
    while (i >= 0) {
        unsigned char c = (unsigned char)buf[i];
        if (c > 13)
            break;
        if (((9217 >> c) & 1) == 0)
            break;
        buf[i--] = 0;
        if (buf == &buf[i + 1])
            break;
    }
    wchar_t *w = mbs_to_wchar(buf);
    char *utf8 = wchar_to_utf8(w);
    lua_pushstring(L, utf8 ? utf8 : "");
    free(utf8);
    free(w);
    return 1;
}

int otogen_lua_is_consonant(lua_State *L)
{
    otogen_ctx *ctx = *(otogen_ctx **)lua_touserdata(L, 1);
    const char *ph = lua_tolstring(L, 2, NULL);
    wchar_t *w = utf8_to_wchar(ph);
    int yes = otogen_is_consonant(ctx->phdict, w);
    free(w);
    lua_settop(L, -3);
    lua_pushboolean(L, yes);
    return 1;
}

wchar_t *otogen_fastlz_to_wchar(const void *src, size_t clen, int dlen)
{
    size_t cap = (size_t)dlen + 512;
    void *buf = calloc(cap, 1);
    memcpy(buf, src, clen);
    size_t n = clen;
    int pass = 2;
    do {
        void *tmp = malloc(cap);
        n = (size_t)fastlz_decompress(buf, (int)n, tmp, (int)cap);
        memcpy(buf, tmp, n);
        free(tmp);
        pass--;
    } while (pass != 0);
    if (strlen((const char *)buf) != n)
        ((char *)buf)[n] = 0;
    int count = MultiByteToWideChar(65001, 0, (LPCCH)buf, -1, NULL, 0);
    wchar_t *wide = (wchar_t *)calloc((size_t)count, sizeof(wchar_t));
    MultiByteToWideChar(65001, 0, (LPCCH)buf, -1, wide, count);
    free(buf);
    return wide;
}

const char *otogen_luac_reader(lua_State *L, void *ud, size_t *sz)
{
    (void)L;
    (void)ud;
    if (g_otogen_luac_fed != 0)
        return NULL;
    g_otogen_luac_fed = 1;
#if MORESAMPLER_ARCH == 64
    *sz = otogen_luac64_size;
    return (const char *)otogen_luac64;
#else
    *sz = otogen_luac_size;
    return (const char *)otogen_luac;
#endif
}

static void push_label_table(lua_State *L, otogen_ctx *ctx, int fi)
{
    otogen_label *lab = ctx->labels[fi];
    wchar_t *base = path_basename(ctx->wav_paths[fi]);

    lua_createtable(L, 0, 0);
    if (lab != NULL) {
        for (int k = 0; k < lab->nseg; k++) {
            lua_createtable(L, 0, 0);
            lua_pushnumber(L, lab->tbegin[k]);
            lua_setfield(L, -2, "tbegin");
            lua_pushnumber(L, lab->tend[k]);
            lua_setfield(L, -2, "tend");
            char *ph = wchar_to_utf8(lab->phoneme[k]);
            lua_pushstring(L, ph ? ph : "");
            free(ph);
            lua_setfield(L, -2, "phoneme");
            if (lab->grapheme[k] != NULL) {
                char *gr = wchar_to_utf8(lab->grapheme[k]);
                lua_pushstring(L, gr);
                free(gr);
            } else {
                lua_pushnil(L);
            }
            lua_setfield(L, -2, "grapheme");
            lua_seti(L, -2, k + 1);
        }
    }
    char *name = wchar_to_utf8(base);
    lua_pushstring(L, name ? name : "");
    lua_setfield(L, -2, "name");
    free(name);
    free(base);
}

static int load_index_csv(otogen_ctx *ctx, const wchar_t *dir, const wchar_t *text)
{
    wchar_t *pref = path_dir_prefix(dir);
    int nline = count_lines(text);
    int n = 0;
    signed int len = (signed int)wcslen(text);
    signed int i = 0;

    ctx->wav_paths = (wchar_t **)calloc((size_t)nline, sizeof(wchar_t *));
    ctx->aliases = (wchar_t **)calloc((size_t)nline, sizeof(wchar_t *));
    while (i < len) {
        signed int start = i;
        signed int end;
        wchar_t *line;
        size_t comma = 0;

        if (text[i] == 13 || text[i] == 10)
            break;
        end = i;
        while (end < len && text[end] != 10 && text[end] != 13)
            end++;
        line = wstr_ndup(&text[start], (size_t)(end - start));
        while (comma < (size_t)(end - start) && line[comma] != L',')
            comma++;
        if (comma >= (size_t)(end - start)) {
            free(line);
            break;
        }
        {
            wchar_t *fn = wstr_ndup(line, comma);
            ctx->wav_paths[n] = wstr_concat(pref, fn);
            ctx->aliases[n] = wstr_dup(line + comma + 1);
            free(fn);
            n++;
        }
        free(line);
        i = end;
        while (i < len && (text[i] == 10 || text[i] == 13))
            i++;
        if (i >= len)
            break;
    }
    ctx->nfile = n;
    free(pref);
    return n;
}

static int greedy_unmatched(otogen_g2p *g, const wchar_t *alias)
{
    signed int n = (signed int)wcslen(alias);
    int pos = 0;

    if (n <= 0)
        return 1;
    while (pos < n) {
        int take = (int)g->max_glen;
        int found = 0;
        if (n - pos < take)
            take = n - pos;
        while (take > 0) {
            wchar_t *s = wstr_ndup(&alias[pos], (size_t)take);
            otogen_hnode *h = otogen_hash_find(g->map, s);
            free(s);
            if (h != NULL) {
                pos += (int)wcslen(h->key);
                found = 1;
                break;
            }
            take--;
        }
        if (!found)
            return 1;
    }
    return 0;
}

int main_oto_generation(int argc, char **argv)
{
    (void)argc;
    lua_State *L = luaL_newstate();
    luaL_openlibs(L);
    luaL_checkversion(L);
    lua_createtable(L, 0, 3);
    luaL_setfuncs(L, otogen_ext_regs, 0);
    lua_setglobal(L, "ext");

    int err = lua_load(L, otogen_luac_reader, NULL, "lmain", "binary");
    if (err != 0)
        goto interp_fail;
    err = lua_pcall(L, 0, 0, 0);
    if (err != 0) {
        lua_settop(L, -2);
        goto interp_fail;
    }

    log_tee_printf("=== Moresampler oto.ini Generation Mode ===\n");
    wchar_t *src_dir = mbs_to_wchar(argv[1]);
    wchar_t *oto_path;
    {
        wchar_t *pref = path_dir_prefix(src_dir);
        oto_path = wstr_concat(pref, L"oto.ini");
        free(pref);
    }

    g_g2p_hiragana = otogen_fastlz_to_wchar(g2p_hiragana_flz, g2p_hiragana_flz_size, g2p_hiragana_dlen);
    g_g2p_romaji = otogen_fastlz_to_wchar(g2p_romaji_flz, g2p_romaji_flz_size, g2p_romaji_dlen);
    g_g2p_english = otogen_fastlz_to_wchar(g2p_english_flz, g2p_english_flz_size, g2p_english_dlen);

    otogen_ctx *ctx = otogen_ctx_create();
    ctx->fft_sec = 0.4f;

    log_tee_printf("1) Japanese (hiragana/katakana)\n2) Japanese (romaji)\n3) English (arpasing)\n");
    log_tee_printf("Please select the input format: ");
    char buf[256];
    fgets(buf, 256, stdin);
    const char *lang;
    const char *variant;
    if (buf[0] == '1') {
        ctx->g2p = otogen_g2p_from_csv(g_g2p_hiragana);
        ctx->phdict = otogen_phdict_from_json(k_phoneme_jp_json);
        variant = "hiragana";
        lang = "japanese";
    } else if (buf[0] == '2') {
        ctx->g2p = otogen_g2p_from_csv(g_g2p_romaji);
        ctx->phdict = otogen_phdict_from_json(k_phoneme_jp_json);
        variant = "romaji";
        lang = "japanese";
    } else {
        ctx->g2p = otogen_g2p_from_csv(g_g2p_english);
        ctx->phdict = otogen_phdict_from_json(k_phoneme_en_json);
        ctx->pad_sec = 0.07f;
        variant = "arpasing";
        lang = "english";
    }
    log_tee_printf("Dictionary loaded. Language: %s, variant: %s\n", lang, variant);

    lua_getglobal(L, "config");
    lua_createtable(L, 0, 0);
    lua_pushstring(L, lang);
    lua_setfield(L, -2, "language");
    lua_pushstring(L, variant);
    lua_setfield(L, -2, "variant");
    lua_setfield(L, -2, "format");

    wchar_t *idx_dir = path_dir_prefix(src_dir);
    wchar_t *idx_path = wstr_concat(idx_dir, L"index.csv");
    wchar_t *idx_text = read_file_acp(idx_path);
    if (idx_text == NULL) {
        free(idx_path);
        free(idx_dir);
        if (otogen_scan_wavs(ctx, src_dir) == 0) {
            log_tee_printf("Error: can't access directory.\n");
            getchar();
            return 1;
        }
        log_tee_printf("%d .wav files detected.\n", ctx->nfile);
        if (ctx->nfile <= 14) {
            log_tee_printf("Insufficient data to initialize training. At least 15 .wav files are required.\n");
            getchar();
            return 1;
        }
        if (otogen_aliases_from_wavs(ctx) == 0) {
            getchar();
            return 1;
        }
    } else {
        load_index_csv(ctx, src_dir, idx_text);
        free(idx_text);
        free(idx_path);
        free(idx_dir);
    }

    {
        wchar_t *pref = path_dir_prefix(src_dir);
        wchar_t *oto_in = wstr_concat(pref, L"oto.ini");
        wchar_t *oto_w = read_file_acp(oto_in);
        free(pref);
        free(oto_in);
        if (oto_w != NULL) {
            char *utf8 = wchar_to_utf8(oto_w);
            lua_getglobal(L, "config");
            lua_pushstring(L, utf8);
            lua_setfield(L, -2, "in_oto");
            lua_settop(L, -2);
            free(utf8);
            free(oto_w);
        }
    }
    free(src_dir);

    lua_getglobal(L, "lconfig");
    err = lua_pcall(L, 0, 0, 0);
    if (err != 0) {
        lua_settop(L, -2);
        log_tee_printf("Error: failed to configure oto generation (error code = %d).\n", err);
        goto bug_fail;
    }

    log_tee_printf("Enter a suffix to be inserted after all entries (or leave it blank): ");
    fgets(buf, 256, stdin);
    {
        wchar_t *suf = mbs_to_wchar(buf);
        signed int sl = (signed int)wcslen(suf);
        if (sl > 0) {
            suf[sl - 1] = 0;
            if (sl != 1) {
                char *utf8 = wchar_to_utf8(suf);
                lua_getglobal(L, "config");
                lua_pushstring(L, utf8);
                lua_setfield(L, -2, "suffix");
                lua_settop(L, -2);
                free(utf8);
                for (int i = 0; i < ctx->nfile; i++) {
                    const wchar_t *al = ctx->aliases[i];
                    signed int aln = (signed int)wcslen(al);
                    if (sl - 1 < aln && wcscmp(&al[aln - (sl - 1)], suf) == 0)
                        ctx->aliases[i][aln - (sl - 1)] = 0;
                }
            }
        }
        free(suf);
    }

    log_tee_printf("Listing files containing unidentified characters:\n");
    {
        int bad = 0;
        for (int i = 0; i < ctx->nfile; i++) {
            wchar_t *base = path_basename(ctx->wav_paths[i]);
            if (greedy_unmatched(ctx->g2p, ctx->aliases[i])) {
                char *m = wchar_to_mbs(base);
                log_tee_printf("%s\n", m);
                free(m);
                bad++;
            }
            free(base);
        }
        log_tee_printf("(%d files in total)\n", bad);
        if (bad != 0)
            log_tee_printf("Removing these files may improve estimation accuracy.\n");
    }

    log_tee_printf("Loading and analyzing audio...\n");
    if (otogen_analyze_audio(ctx) == 0)
        goto wait_fail;
    if (ctx->observset == NULL)
        goto wait_fail;

    log_tee_printf("Initializing segmentation...\n");
    otogen_init_segmentation(ctx);
    if (ctx->segset == NULL)
        goto wait_fail;

    log_tee_printf("Initializing model...\n");
    otogen_init_model(ctx, 0);
    ctx->boundary_model = otogen_bmodel_from_table();
    if (ctx->model == NULL)
        goto wait_fail;

    lrh_output_vfloor = 30.0f;
    lrh_duration_vfloor = 0.5f;
    ctx->model->streams[0]->weight = 0.5f;

    log_tee_printf("Training (1st stage)...\n");
    {
        int it = 0;
        FP_TYPE prev = 0.0f;

        for (;;) {
            ++it;
            printf("Iteration %d/%d: ", it, 5);
            fflush(stdout);
            lrh_model_stat *st = lrh_model_stat_from_model(ctx->model);
            FP_TYPE sum = 0.0f;
            if (ctx->nfile > 0) {
                for (int fi = 0; fi < ctx->nfile; fi++) {
                    lrh_observ *ob = ctx->observset->samples[fi];
                    lrh_seg *seg = ctx->segset->samples[fi];

                    seg->time[seg->nseg - 1] = ob->nt;
                    lrh_inference_stprune = (FP_TYPE)ctx->infer_nt[fi];
                    FP_TYPE ll = lrh_estimate_geometric(st, ctx->model, ob, seg);
                    if (ll < -1000000000.0)
                        fprintf(stderr,
                            "Warning: unsuccessful inference on sample %d/%d.\n",
                            fi + 1, ctx->nfile);
                    sum += ll;
                }
            }
            FP_TYPE mean = sum / (FP_TYPE)ctx->nfile;
            printf("Mean log likelihood = %f\n", (double)mean);
            lrh_model_update(ctx->model, st, 1);
            otogen_pool_vowel_gmms(ctx);
            lrh_delete_model_stat(st);
            if ((unsigned)(it - 1) > 1 && prev + 1.0 >= (double)mean) {
                printf("Convergence reached at iteration %d.\n", it);
                break;
            }
            prev = mean;
            if (it == 5)
                break;
        }
    }

    log_tee_printf("Training (2st stage)...\n");
    otogen_train(ctx, 7, 1);

    log_tee_printf("Editing model...\n");
    otogen_edit_silence_floor(ctx);

    lrh_output_vfloor = 1.0f;
    lrh_duration_vfloor = 0.01f;
    puts("Training (3rd stage)...");
    otogen_train(ctx, 6, 1);
    puts("Aligning...");
    otogen_align(ctx, 0);
    puts("Refining alignment...");
    otogen_refine_alignment(ctx, 1);
    otogen_fill_label_times(ctx);

    log_tee_printf("Generating oto.ini...\n");
    lua_getglobal(L, "lmain");
    lua_createtable(L, 0, 0);
    lua_createtable(L, 0, 0);
    lua_pushstring(L, lang);
    lua_setfield(L, -2, "language");
    lua_pushstring(L, variant);
    lua_setfield(L, -2, "variant");
    lua_setfield(L, -2, "format");
    lua_createtable(L, 0, 0);
    for (int i = 0; i < ctx->nfile; i++) {
        push_label_table(L, ctx, i);
        lua_seti(L, -2, i + 1);
    }
    lua_setfield(L, -2, "label_list");
    {
        otogen_ctx **box = (otogen_ctx **)lua_newuserdata(L, sizeof(otogen_ctx *));
        *box = ctx;
    }
    lua_setfield(L, -2, "ptr");
    err = lua_pcall(L, 1, 1, 0);
    if (err != 0) {
        lua_settop(L, -2);
        log_tee_printf("Error: failed to generate oto (error code = %d).\n", err);
        goto bug_fail;
    }

    const char *result = lua_tolstring(L, -1, NULL);
    wchar_t *result_w = utf8_to_wchar(result);
    FILE *out = _wfopen(oto_path, L"w");
    if (out == NULL) {
        log_tee_printf("Error: can't write oto.ini.\n");
        goto wait_fail;
    }
    log_tee_printf("Saving...\n");
    char *result_mbs = wchar_to_mbs(result_w);
    fputs(result_mbs ? result_mbs : "", out);
    fclose(out);
    free(result_mbs);
    free(result_w);
    free(oto_path);
    system("pause");
    otogen_ctx_free(ctx);
    free(g_g2p_hiragana);
    lua_close(L);
    return 0;

interp_fail:
    log_tee_printf("Error: failed to initialize interpreter (error code = %d).\n", err);
bug_fail:
    log_tee_printf("It's likely that you ran into a bug. Please send a bug report to the author <k.hua.kanru [at] ieee.org>\n");
    log_tee_printf("Oto generation terminated.\n");
    return 1;
wait_fail:
    getchar();
    return 1;
}
