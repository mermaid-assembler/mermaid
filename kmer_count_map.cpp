#include "config.h"
#include "kmer.h"
#include "kmer_count_map.h"
#include "utils.h"
#include "contig_store.h"
#include "kmer_ext_map.h"

#define FILTER_ON 0

KmerCountMap::KmerCountMap(k_t k)
    : k(k),
    sb_filter(INITIAL_CAPACITY, 0.001, 0,
            (bloom_filter_hash_func_t) kmer_hash_K),
    hash_map(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k))
{
}

void KmerCountMap::insert(qekmer_t* qekmer)
{
#if FILTER_ON
    if (!sb_filter.check(qekmer->kmer)) {
        sb_filter.add(qekmer->kmer);
    } else {
#endif
        hash_map.try_insert(qekmer->kmer);
        if (qekmer->lqual > Config::Q_MIN)
            hash_map.map[qekmer->kmer].lquals[qekmer->exts.left]++;
        if (qekmer->rqual > Config::Q_MIN)
            hash_map.map[qekmer->kmer].rquals[qekmer->exts.right]++;
#if FILTER_ON
    }
#endif
}

void KmerCountMap::trim(KmerExtMap& kmer_ext_map)
{
    for (map_type_t::iterator it = hash_map.map.begin();
            it != hash_map.map.end();
            it++) {
        kmer_t kmer = it->first;
        ext_map_t ext_map = it->second.ext_map(Config::D_MIN);

        if (ext_map.valid())
            kmer_ext_map.insert(kmer, ext_map);
        else
            free(kmer);
    }
}
