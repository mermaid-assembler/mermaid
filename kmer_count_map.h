#ifndef _KMER_COUNT_MAP_H_
#define _KMER_COUNT_MAP_H_

#include <vector>

#include "kmer.h"
#include "hash_map.h"
#include "scalable_bloom_filter.h"
#include "contig.h"
#include "contig_store.h"

class KmerExtMap;

class KmerCountMap {
public:
    KmerCountMap(k_t k);

    /* Insert qekmer into the kmer hash map. */
    void insert(qekmer_t* qekmer);
    
    /* Trim kmers which don't appear at least D_MIN times. */
    void trim(KmerExtMap& kmer_ext_map);

protected:
    typedef HashMap<kmer_t, qual_counts_t>::map_type_t map_type_t;

    k_t k;
    ScalableBloomFilter sb_filter;
    HashMap<kmer_t, qual_counts_t> hash_map;

    /* FIXME - Change this initial capacity using preprocessing step. */
    static const int INITIAL_CAPACITY = 100000;
};

#endif /* _KMER_COUNT_MAP_H_ */
