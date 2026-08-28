#include "otogen.h"

#include "../logging.h"
#include "../utils.h"

static otogen_plist *plist_new(const wchar_t *ph)
{
    otogen_plist *n = (otogen_plist *)malloc(sizeof(otogen_plist));
    n->phoneme = wstr_dup(ph);
    n->unused = 0;
    n->next = NULL;
    return n;
}

static void plist_free(otogen_plist *n)
{
    while (n != NULL) {
        otogen_plist *nx = n->next;
        free(n->phoneme);
        free(n);
        n = nx;
    }
}

static otogen_plist *parse_phoneme_chain(const wchar_t *s)
{
    const wchar_t *p = s;
    otogen_plist *head = NULL;
    otogen_plist **tail = &head;

    while (*p) {
        const wchar_t *start = p;
        while (*p && *p != L',')
            p++;
        {
            wchar_t *tok = wstr_ndup(start, (size_t)(p - start));
            *tail = plist_new(tok);
            free(tok);
        }
        tail = &(*tail)->next;
        if (*p == L',')
            p++;
    }
    return head;
}

otogen_g2p *otogen_g2p_from_csv(const wchar_t *text)
{
    signed int n = (signed int)wcslen(text);
    int nline;
    if (n <= 0) {
        nline = 0;
    } else {
        int breaks = 0;
        signed int p = 0;
        wchar_t c = *text;
        while (1) {
            while (c != 10 && c != 13) {
                if (n <= ++p)
                    goto counted;
                c = text[p];
            }
            if (n <= p)
                break;
            while (c == 10 || c == 13) {
                if (n == ++p)
                    goto counted;
                c = text[p];
            }
            ++breaks;
        }
    counted:
        nline = breaks + 1;
    }

    otogen_g2p *g = (otogen_g2p *)malloc(sizeof(otogen_g2p));
    g->map = otogen_hash_create(nline > 0 ? nline : 1);
    g->chains = (otogen_plist **)malloc((size_t)nline * sizeof(otogen_plist *));
    g->max_glen = 0;
    int line_i = 0;
    signed int i = 0;
    while (i < n) {
        signed int start = i;
        if (text[i] == 13 || text[i] == 10) {
            while (i < n && (text[i] == 10 || text[i] == 13))
                i++;
            continue;
        }
        while (i < n && text[i] != 10 && text[i] != 13)
            i++;
        signed int len = i - start;
        wchar_t *line = wstr_ndup(&text[start], (size_t)len);
        wchar_t *comma = wcschr(line, L',');
        wchar_t *grapheme;
        otogen_plist *chain;
        if (comma == NULL) {
            grapheme = wstr_dup(line);
            chain = NULL;
        } else {
            grapheme = wstr_ndup(line, (size_t)(comma - line));
            chain = parse_phoneme_chain(comma + 1);
        }
        otogen_hash_insert(g->map, grapheme, line_i);
        g->chains[line_i] = chain;
        size_t glen = wcslen(grapheme);
        if (g->max_glen < glen)
            g->max_glen = glen;
        free(grapheme);
        free(line);
        line_i++;
        while (i < n && (text[i] == 10 || text[i] == 13))
            i++;
    }
    return g;
}

void otogen_g2p_free(otogen_g2p *g)
{
    if (g == NULL)
        return;
    if (g->chains != NULL && g->map != NULL) {
        for (int i = 0; i < g->map->count; i++)
            plist_free(g->chains[i]);
        free(g->chains);
    }
    otogen_hash_free(g->map);
    free(g);
}

void otogen_label_free(otogen_label *lab)
{
    if (lab == NULL)
        return;
    for (int i = 0; i < lab->nseg; i++) {
        free(lab->phoneme[i]);
        free(lab->grapheme[i]);
    }
    free(lab->tbegin);
    free(lab->tend);
    free(lab->phoneme);
    free(lab->grapheme);
    free(lab);
}

otogen_label *otogen_g2p_apply(otogen_g2p *g, const wchar_t *text)
{
    signed int n = (signed int)wcslen(text);
    otogen_plist sil1;
    sil1.phoneme = wstr_dup(L"sil1");
    sil1.unused = -1;
    sil1.next = NULL;
    otogen_plist *head = &sil1;
    otogen_plist *cur = head;
    int pos = 0;
    while (pos < n) {
        int take = (int)g->max_glen;
        otogen_hnode *hit = NULL;

        if (n - pos < take)
            take = n - pos;
        while (take > 0) {
            wchar_t *slice = wstr_ndup(&text[pos], (size_t)take);
            hit = otogen_hash_find(g->map, slice);
            free(slice);
            if (hit != NULL)
                break;
            take--;
        }
        if (hit == NULL) {
            wchar_t *one = wstr_ndup(&text[pos], 1);
            char *mbs = wchar_to_mbs(one);
            fprintf(stderr, "Unidentified grapheme: %s (skipped)\n", mbs);
            free(mbs);
            free(one);
            pos++;
            continue;
        }
        {
            otogen_plist *nnode = plist_new(hit->key);
            nnode->unused = hit->index;
            cur->next = nnode;
            cur = nnode;
            pos += (int)wcslen(hit->key);
        }
    }
    otogen_plist sil2;
    sil2.phoneme = wstr_dup(L"sil2");
    sil2.unused = -2;
    sil2.next = NULL;
    cur->next = &sil2;

    int nseg = 0;
    otogen_plist *walk = head;
    while (walk != NULL) {
        if (walk->unused < 0)
            nseg++;
        else {
            otogen_plist *ch = g->chains[walk->unused];
            while (ch != NULL) {
                nseg++;
                ch = ch->next;
            }
        }
        walk = walk->next;
    }

    otogen_label *lab = (otogen_label *)malloc(sizeof(otogen_label));
    lab->tbegin = (FP_TYPE *)calloc((size_t)nseg, sizeof(FP_TYPE));
    lab->tend = (FP_TYPE *)calloc((size_t)nseg, sizeof(FP_TYPE));
    lab->phoneme = (wchar_t **)calloc((size_t)nseg, sizeof(wchar_t *));
    lab->grapheme = (wchar_t **)calloc((size_t)nseg, sizeof(wchar_t *));
    lab->nseg = nseg;
    int out_i = 0;
    walk = head;
    while (walk != NULL) {
        if (walk->unused < 0) {
            lab->phoneme[out_i] = wstr_dup(walk->phoneme);
            lab->grapheme[out_i] = NULL;
            out_i++;
        } else {
            otogen_plist *ch = g->chains[walk->unused];
            if (ch != NULL)
                lab->grapheme[out_i] = wstr_dup(walk->phoneme);
            while (ch != NULL) {
                lab->phoneme[out_i] = wstr_dup(ch->phoneme);
                out_i++;
                ch = ch->next;
            }
        }
        walk = walk->next;
    }

    free(sil1.phoneme);
    free(sil2.phoneme);
    walk = head->next;
    while (walk != NULL && walk != &sil2) {
        otogen_plist *nx = walk->next;
        free(walk->phoneme);
        free(walk);
        walk = nx;
    }
    return lab;
}

otogen_phdict *otogen_phdict_from_json(const char *json)
{
    cJSON *root = cJSON_Parse(json);
    if (root == NULL)
        return NULL;
    cJSON *phoneme = cJSON_GetObjectItem(root, "phoneme-list");
    cJSON *cons = cJSON_GetObjectItem(root, "consonant-list");
    cJSON *vowel = cJSON_GetObjectItem(root, "vowel-list");
    cJSON *voiced = cJSON_GetObjectItem(root, "voiced-list");
    cJSON *unvoiced = cJSON_GetObjectItem(root, "unvoiced-list");
    cJSON *semi = cJSON_GetObjectItem(root, "semi-list");
    cJSON *silence = cJSON_GetObjectItem(root, "silence-list");
    if (phoneme == NULL || cons == NULL || vowel == NULL || voiced == NULL
        || unvoiced == NULL || semi == NULL || silence == NULL) {
        cJSON_Delete(root);
        return NULL;
    }
    int n = cJSON_GetArraySize(phoneme);
    otogen_phdict *d = (otogen_phdict *)malloc(sizeof(otogen_phdict));
    d->names = (wchar_t **)malloc((size_t)n * sizeof(wchar_t *));
    d->all = otogen_hash_create(n > 0 ? n : 1);
    int i = 0;
    for (cJSON *it = phoneme->child; it != NULL; it = it->next) {
        wchar_t *w = mbs_to_wchar(it->valuestring);
        otogen_hash_insert(d->all, w, i);
        d->names[i++] = w;
    }
    d->consonant = otogen_hash_from_cjson_array(cons, d->all);
    d->vowel = otogen_hash_from_cjson_array(vowel, d->all);
    d->voiced = otogen_hash_from_cjson_array(voiced, d->all);
    d->unvoiced = otogen_hash_from_cjson_array(unvoiced, d->all);
    d->semi = otogen_hash_from_cjson_array(semi, d->all);
    d->silence = otogen_hash_from_cjson_array(silence, d->all);
    cJSON_Delete(root);
    return d;
}

void otogen_phdict_free(otogen_phdict *d)
{
    if (d == NULL)
        return;
    if (d->names != NULL && d->all != NULL) {
        for (int i = 0; i < d->all->count; i++)
            free(d->names[i]);
        free(d->names);
    }
    otogen_hash_free(d->all);
    otogen_hash_free(d->consonant);
    otogen_hash_free(d->vowel);
    otogen_hash_free(d->voiced);
    otogen_hash_free(d->unvoiced);
    otogen_hash_free(d->semi);
    otogen_hash_free(d->silence);
    free(d);
}

int otogen_is_consonant(otogen_phdict *d, const wchar_t *ph)
{
    return d != NULL && otogen_hash_has(d->consonant, ph);
}

int otogen_is_vowel(otogen_phdict *d, const wchar_t *ph)
{
    return d != NULL && otogen_hash_has(d->vowel, ph);
}

int otogen_is_voiced(otogen_phdict *d, const wchar_t *ph)
{
    return d != NULL && otogen_hash_has(d->voiced, ph);
}

int otogen_is_unvoiced(otogen_phdict *d, const wchar_t *ph)
{
    return d != NULL && otogen_hash_has(d->unvoiced, ph);
}

int otogen_is_semi(otogen_phdict *d, const wchar_t *ph)
{
    return d != NULL && otogen_hash_has(d->semi, ph);
}

int otogen_is_silence(otogen_phdict *d, const wchar_t *ph)
{
    return d != NULL && otogen_hash_has(d->silence, ph);
}
