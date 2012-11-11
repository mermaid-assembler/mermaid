#include <cassert>
#include <cstdio>

#include "kmer.h"
#include "config.h"
#include "utils.h"

int cmp_kmer(kmer_t kmer1, kmer_t kmer2, k_t len, k_t start1, k_t start2)
{
    base b1, b2;
    for_base_in_kmer_from(b1, kmer1, len, start1) {
        b2 = get_base(kmer2, start2 + b1_i_);
        if (b1 < b2) return -1;
        else if (b1 > b2) return 1;
    } end_for;

    return 0;
}

void str2kmer(kmer_t kmer, const char* str, k_t len)
{
    for (k_t i = 0; i < len; i++) {
        set_base(kmer, i, char2base(str[i]));
    }
}

void kmer2str(char* str, kmer_t kmer, k_t len, k_t start)
{
    base b;
    for_base_in_kmer_from(b, kmer, len, start) {
        str[b_i_] = base2char(b);
    } end_for;
    str[len] = '\0';
}

void revcmp_kmer(kmer_t dst, kmer_t src, k_t len)
{
    base b;
    for_base_in_kmer(b, src, len) {
        set_base(dst, (len - 1) - b_i_, inv_base(b));
    } end_for;
}

void fprint_kmer(FILE* f, kmer_t kmer, k_t len, k_t start)
{
    char str[len + 1];
    kmer2str(str, kmer, len, start);
    str[len] = '\0';

    fprintf(f, "%s", str);
}

void fprintln_kmer(FILE* f, kmer_t kmer, k_t len, k_t start)
{
    fprint_kmer(f, kmer, len, start);
    fprintf(f, "\n");
}

void println_kmer(kmer_t kmer, k_t len, k_t start)
{
    fprintln_kmer(stdout, kmer, len, start);
}

bool kmer_matches_str(kmer_t kmer, const char* str, k_t len, k_t from)
{
    kmer_a to_match[kmer_size(len)];
    str2kmer(to_match, str, len);
    return cmp_kmer(kmer, to_match, len, from) == 0;
}


bool kmer_matches_str_or_revcmp(kmer_t kmer, const char* str, k_t len)
{
    kmer_a revcmp[kmer_size(len)];
    revcmp_kmer(revcmp, kmer, len);
    return kmer_matches_str(kmer, str, len) || cmp_kmer(revcmp, kmer, len) == 0;
}

void hashlittle2(const void *key, size_t length, uint32_t *pc, uint32_t *pb);

size_t kmer_hash(size_t seed, kmer_t kmer, k_t len)
{
    size_t ret_val;
    size_t byte_len = len / BASES_PER_BYTE;
    uint32_t pc = (uint32_t) (seed >> 32);
    uint32_t pb = (uint32_t) seed;

    hashlittle2(kmer, byte_len, &pc, &pb);
    ret_val = pc + (((size_t) pb) << 32);

    base b;
    for_base_in_kmer_from(b, kmer, len % BASES_PER_BYTE, byte_len * BASES_PER_BYTE) {
        /* This is from djb2. */
        ret_val = ((ret_val << 5) + ret_val) + (size_t) b;
    } end_for;

    return ret_val;
}

size_t kmer_hash_K(size_t seed, kmer_t kmer)
{
    return kmer_hash(seed, kmer, Config::K);
}

size_t kmer_hash_simple(kmer_t kmer)
{
    return kmer_hash(0, kmer, Config::K);
}

bool kmer_eq(kmer_t x, kmer_t y, k_t len)
{
    return cmp_kmer(x, y, len) == 0;
}

bool kmer_eq_K(kmer_t x, kmer_t y)
{
    return cmp_kmer(x, y, Config::K) == 0;
}

bool is_canonical_kmer(kmer_t kmer, k_t len)
{
    kmer_a revcmp[kmer_size(len)];
    revcmp_kmer(revcmp, kmer, len);
    return cmp_kmer(kmer, revcmp, len) <= 0;
}

bool validate_kmer(kmer_t kmer, k_t len)
{
    base b;
    for_base_in_kmer(b, kmer, len) {
        if (!BASE::valid_base(b))
            return false;
    } end_for;
    return true;
}

void assert_kmer(kmer_t kmer, k_t len)
{
    if (!validate_kmer(kmer, len))
        panic("%s\n", __func__);
}
