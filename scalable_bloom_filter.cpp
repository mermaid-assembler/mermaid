#include <cstdlib>
#include <cmath>
#include <vector>

#include "scalable_bloom_filter.h"

using namespace std;

static const size_t SCALE = 2;
static const float RATIO = 0.9;

ScalableBloomFilter::ScalableBloomFilter(size_t initial_capacity,
        float error_rate, size_t seed, bloom_filter_hash_func_t hash)
    : filters(), initial_capacity(initial_capacity), error_rate(error_rate),
    scale(SCALE), ratio(RATIO), seed(seed), hash(hash)
{
    filters.push_back(new BloomFilter(initial_capacity, error_rate, seed, hash));
}

ScalableBloomFilter::~ScalableBloomFilter()
{
    while (!filters.empty()) {
        delete filters.back();
        filters.pop_back();
    }
}

void ScalableBloomFilter::add(void* key)
{
    BloomFilter* filter = filters.back();
    if (filter->full()) {
        filter = new BloomFilter(initial_capacity * pow(scale, (float) filters.size()),
                error_rate * pow(ratio, (float) filters.size()),
                filters.size() + seed, hash);
        filters.push_back(filter);
    }
    filter->add(key);
}

bool ScalableBloomFilter::check(void* key)
{
    for (vector<BloomFilter*>::iterator it = filters.begin();
            it < filters.end();
            it++) {
        if ((*it)->check(key))
            return true;
    }
    return false;
}
