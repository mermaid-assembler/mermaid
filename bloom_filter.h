#ifndef _BLOOM_FILTER_H_
#define _BLOOM_FILTER_H_

#include <boost/cstdint.hpp>

typedef size_t (*bloom_filter_hash_func_t)(size_t seed, void* key);

class BloomFilter {
public:
    BloomFilter(size_t capacity, float error_rate, size_t seed, bloom_filter_hash_func_t hash);
    ~BloomFilter();

    /* Adds the key to the filter. */
    void add(void* key);
    /* Returns true if key was found in the filter; false otherwise. */
    bool check(void* key);

    bool full(void) { return count >= capacity; }

protected:
    uint8_t* array;
    size_t capacity;
    size_t nfuncs;
    size_t slice_size;          /* Size of each partition. */
    size_t count;
    size_t seed;
    bloom_filter_hash_func_t hash;
};


#endif /* _BLOOM_FILTER_H_ */
