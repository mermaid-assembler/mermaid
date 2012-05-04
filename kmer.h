#ifndef _KMER_H_
#define _KMER_H_

#include <cstdio>

#include <boost/cstdint.hpp>

#include "config.h"
#include "base.h"
#include "utils.h"

using BASE::base;
using BASE::char2base;
using BASE::base2char;
using BASE::inv_base;
using BASE::valid_base;

/* This must be a power of 2. */
const uint8_t BITS_PER_BASE = 2;
const uint8_t BASES_PER_BYTE = 8 / BITS_PER_BASE;

typedef uint8_t kmer_base_t;
typedef uint8_t* kmer_t;
typedef uint8_t kmer_a;     /* This type should be used for declaring arrays
                               that will hold kmers. */

typedef struct {
    uint8_t left : 4;
    uint8_t right : 4;

    base left_ext(void)
    {
        if (!is_uu())
            panic("ext_map is not uu\n");

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (left & 1 << i)
                return (base) i;
        }

        panic("ext_map is not uu\n");
    }

    base right_ext(void)
    {
        if (!is_uu())
            panic("ext_map is not uu\n");

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (right & 1 << i)
                return (base) i;
        }

        /* To shut gcc up.... */
        panic("ext_map is not uu\n");
    }

    void invalidate(void)
    {
        left = 0;
        right = 0;
    }

    bool valid(void)
    {
        return left || right;
    }

    bool is_uu(void)
    {
        uint8_t valid_left_bases = 0;
        uint8_t valid_right_bases = 0;

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (left & 1 << i)
                valid_left_bases++;
            if (right & 1 << i)
                valid_right_bases++;
        }

        return valid_left_bases == 1 && valid_right_bases == 1;
    }
} ext_map_t;

typedef struct {
    count_t lquals[BASE::NUM_BASES];
    count_t rquals[BASE::NUM_BASES];

    ext_map_t ext_map(uint8_t d_min)
    {
        ext_map_t ext_map = {0, 0};
        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (lquals[i] >= d_min)
                ext_map.left |= 1 << i;
            if (rquals[i] >= d_min)
                ext_map.right |= 1 << i;
        }
        return ext_map;
    }
} qual_counts_t;

typedef struct {
    uint8_t left : BITS_PER_BASE;
    uint8_t right : BITS_PER_BASE;
} exts_t;

typedef struct {
    qual_t lqual;
    qual_t rqual;
    exts_t exts;
    kmer_a kmer[0];
} __attribute__((packed)) qekmer_t;

#define kmer_size(len) \
    (((len) + BASES_PER_BYTE - 1) / BASES_PER_BYTE)

#define qekmer_size(len) \
    (sizeof(qekmer_t) + kmer_size(len))

/* WARNING: This function does not check its parameters for correct input. */
#define get_base(kmer, i) \
    ((base) (((1 << BITS_PER_BASE) - 1) & \
        ((kmer)[(i) / BASES_PER_BYTE] >> (BITS_PER_BASE * ((i) % BASES_PER_BYTE)))))

/* WARNING: This function does not check its parameters for correct input. */
#define set_base(kmer, i, b)                                            \
    do {                                                                    \
        kmer_a mask = (1 << BITS_PER_BASE) - 1;                            \
        uint8_t shift = BITS_PER_BASE * ((i) % BASES_PER_BYTE);           \
        (kmer)[(i) / BASES_PER_BYTE] &= ~(mask << shift);                   \
        (kmer)[(i) / BASES_PER_BYTE] |= (mask & (b)) << shift;                \
    } while (0)


#define for_base_in_kmer(b, kmer, len)                                \
    do {                                                              \
        for (k_t b ## _i_ = 0; b ## _i_ < (len); b ## _i_++) {     \
            b = get_base((kmer), b ## _i_);

#define for_base_in_kmer_from(b, kmer, len, from)                   \
    do {                                                            \
        for (k_t b ## _i_ = 0; b ## _i_ < (len); b ## _i_++) {   \
            b = get_base((kmer), (from) + b ## _i_);

#define end_for } } while (0)

int cmp_kmer(kmer_t kmer1, kmer_t kmer2, k_t len, k_t start1 = 0, k_t start2 = 0);

void str2kmer(kmer_t kmer, const char* str, k_t len);
void kmer2str(char* str, kmer_t kmer, k_t len, k_t start = 0);

void revcmp_kmer(kmer_t dst, kmer_t src, k_t len);

void fprint_kmer(FILE* f, kmer_t kmer, k_t len, k_t start = 0);
void fprintln_kmer(FILE* f, kmer_t kmer, k_t len, k_t start = 0);
void println_kmer(kmer_t kmer, k_t len, k_t start = 0);

size_t kmer_hash(size_t seed, kmer_t kmer, k_t len);
/* Version of kmer_hash in which len is assumed to be K from config.h */
size_t kmer_hash_K(size_t seed, kmer_t kmer);
/* Version of kmer_hash in which len is assumed to be K from config.h
 * and seed = 0 */
size_t kmer_hash_simple(kmer_t kmer);

bool kmer_eq(kmer_t x, kmer_t y, k_t len);
/* Version of kmer_eq in which len is assumed to be K from config.h */
bool kmer_eq_K(kmer_t x, kmer_t y);

bool is_canonical_kmer(kmer_t kmer, k_t len);

/*****************************************************************************
 * The following are useful debugging functions. NOT meant for production use.
 *****************************************************************************/
bool kmer_matches_str(kmer_t kmer, const char* str, k_t len, k_t from = 0);
bool kmer_matches_str_or_revcmp(kmer_t kmer, const char* str, k_t len);

bool validate_kmer(kmer_t kmer, k_t len);

/* Validates kmer and on failure raises an error. */
void assert_kmer(kmer_t kmer, k_t len);

#endif /* _KMER_H_ */
