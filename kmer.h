#ifndef _KMER_H_
#define _KMER_H_

#include <cstdio>

#include <boost/cstdint.hpp>

#include "config.h"
#include "base.h"

using BASE::base;
using BASE::char2base;
using BASE::base2char;
using BASE::inv_base;

/* This must be a power of 2. */
const uint8_t BITS_PER_BASE = 4;
const uint8_t BASES_PER_BYTE = 8 / BITS_PER_BASE;

typedef uint8_t kmer_base_t;
typedef uint8_t* kmer_t;
typedef uint8_t kmer_a;     /* This type should be used for declaring arrays
                               that will hold kmers. */

typedef uint8_t qual_t;
typedef uint32_t count_t;

typedef struct {
    qual_t lqual;
    qual_t rqual;
    struct {
        uint8_t left : 4;
        uint8_t right : 4;
    } exts;
    kmer_a kmer[0];
} __attribute__((packed)) qekmer_t;

/* Packed structure for representing UFX for both left and right */
typedef struct extensions {
    union
    {
        struct
        {
            unsigned int lT:1;
            unsigned int lG:1;
            unsigned int lC:1;
            unsigned int lA:1;
            unsigned int rT:1;
            unsigned int rG:1;
            unsigned int rC:1;
            unsigned int rA:1;
        };
        uint8_t ext;
    };
} extensions_t;

/* Combined extensions_t and kmer */
typedef struct ekmer {
    extensions_t ext;
    kmer_a kmer[0];
} __attribute__((packed)) ekmer_t;

#define kmer_size(len) \
    (((len) + BASES_PER_BYTE - 1) / BASES_PER_BYTE)

#define qekmer_size(len) \
    (sizeof(qekmer_t) + kmer_size(len))

#define ekmer_size(len) \
    (sizeof(ekmer_t) + kmer_size(len))

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

/* The following are useful debugging functions. NOT meant for production use. */
bool kmer_matches_str(kmer_t kmer, const char* str, k_t len, k_t from = 0);

size_t kmer_hash(size_t seed, kmer_t kmer, k_t len);
/* Version of kmer_hash in which len is assumed to be K from config.h */
size_t kmer_hash_K(size_t seed, kmer_t kmer);
/* Version of kmer_hash in which len is assumed to be K from config.h
 * and seed = 0 */
size_t kmer_hash_simple(kmer_t kmer);

bool kmer_eq(kmer_t x, kmer_t y, k_t len);
/* Version of kmer_eq in which len is assumed to be K from config.h */
bool kmer_eq_K(kmer_t x, kmer_t y);

bool validate_kmer(kmer_t kmer, k_t len);

#endif /* _KMER_H_ */
