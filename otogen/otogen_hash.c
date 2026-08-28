#include "otogen.h"
#include "../utils.h"

unsigned otogen_hash_wchar(const wchar_t *s, int nbucket)
{
    signed int n = (signed int)wcslen(s);
    unsigned h = 0;

    if (n <= 0 || nbucket <= 0)
        return 0;
    for (signed int i = 0; i < n; i++)
        h = (unsigned)s[i] ^ (33 * h);
    return h % (unsigned)nbucket;
}

otogen_hash *otogen_hash_create(int nbucket)
{
    otogen_hash *h = (otogen_hash *)malloc(sizeof(otogen_hash));
    h->count = 0;
    h->nbucket = nbucket > 0 ? nbucket : 1;
    h->head = (otogen_hnode **)calloc((size_t)h->nbucket, sizeof(otogen_hnode *));
    h->tail = (otogen_hnode **)calloc((size_t)h->nbucket, sizeof(otogen_hnode *));
    return h;
}

void otogen_hash_insert(otogen_hash *h, const wchar_t *key, int index)
{
    unsigned b = otogen_hash_wchar(key, h->nbucket);
    otogen_hnode *n = (otogen_hnode *)malloc(sizeof(otogen_hnode));
    n->key = wstr_dup(key);
    n->index = index;
    n->next = NULL;
    if (h->head[b] != NULL) {
        h->tail[b]->next = n;
        h->tail[b] = n;
    } else {
        h->head[b] = n;
        h->tail[b] = n;
    }
    h->count++;
}

otogen_hnode *otogen_hash_find(otogen_hash *h, const wchar_t *key)
{
    if (h == NULL || key == NULL)
        return NULL;
    unsigned b = otogen_hash_wchar(key, h->nbucket);
    otogen_hnode *n = h->head[b];
    while (n != NULL) {
        if (wcscmp(n->key, key) == 0)
            return n;
        n = n->next;
    }
    return NULL;
}

int otogen_hash_has(otogen_hash *h, const wchar_t *key)
{
    return otogen_hash_find(h, key) != NULL;
}

void otogen_hash_free(otogen_hash *h)
{
    if (h == NULL)
        return;
    for (int i = 0; i < h->nbucket; i++) {
        otogen_hnode *n = h->head[i];
        while (n != NULL) {
            otogen_hnode *nx = n->next;
            free(n->key);
            free(n);
            n = nx;
        }
    }
    free(h->head);
    free(h->tail);
    free(h);
}

otogen_hash *otogen_hash_from_cjson_array(cJSON *arr, otogen_hash *all)
{
    int n = cJSON_GetArraySize(arr);
    otogen_hash *h = otogen_hash_create(n > 0 ? n : 1);

    for (cJSON *it = arr->child; it != NULL; it = it->next) {
        wchar_t *w = mbs_to_wchar(it->valuestring);
        otogen_hnode *src = all != NULL ? otogen_hash_find(all, w) : NULL;
        otogen_hash_insert(h, w, src != NULL ? src->index : h->count);
        free(w);
    }
    return h;
}
