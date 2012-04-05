#include <cstdlib>
#include <climits>
#include <cmath>

#include <boost/cstdint.hpp>

#include "bloom_filter.h"

#define set_bit(a, i) ((a)[(i)/CHAR_BIT] |= (1 << ((i) % CHAR_BIT)))
#define get_bit(a, i) ((a)[(i)/CHAR_BIT] & (1 << ((i) % CHAR_BIT)))

BloomFilter::BloomFilter(size_t capacity, float error_rate, size_t seed, bloom_filter_hash_func_t hash)
    : capacity(capacity), count(0), seed(seed), hash(hash)
{
    nfuncs = (size_t) ceil(log2(1 / error_rate));
    slice_size = (size_t) ceil(capacity * log(2));
    /* FIXME - Using CHAR_BIT may tie us to using only gcc. */
    array = (uint8_t*) calloc((nfuncs * slice_size + CHAR_BIT - 1) / CHAR_BIT, sizeof(uint8_t));
}

BloomFilter::~BloomFilter()
{
    free(array);
}

void BloomFilter::add(void* key)
{
    for (size_t i = 0; i < nfuncs; i++) {
        set_bit(array, hash(i + seed, key) % slice_size + i * slice_size);
    }
    count++;
}

bool BloomFilter::check(void* key)
{
    for (size_t i = 0; i < nfuncs; i++) {
        if (!get_bit(array, hash(i + seed, key) % slice_size + i * slice_size))
            return false;
    }
    return true;
}
