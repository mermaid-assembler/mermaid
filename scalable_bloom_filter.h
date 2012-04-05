#ifndef _SCALABLE_BLOOM_FILTER_H_
#define _SCALABLE_BLOOM_FILTER_H_

#include <vector>

#include "bloom_filter.h"

class ScalableBloomFilter {
public:
    ScalableBloomFilter(size_t initial_capacity, float error_rate, size_t seed, bloom_filter_hash_func_t hash);
    ~ScalableBloomFilter();

    /* Adds the key to the filter. */
    void add(void* key);
    /* Returns true if key was found in the filter; false otherwise. */
    bool check(void* key);

protected:
    std::vector<BloomFilter*> filters;
    size_t initial_capacity;
    float error_rate;
    float scale;
    float ratio;
    size_t seed;
    bloom_filter_hash_func_t hash;
};


#endif /* _SCALABLE_BLOOM_FILTER_H_ */
