#include "config.h"
#include "kmer.h"
#include "kmer_store.h"
#include "utils.h"
#include "contig_store.h"
#include "trimmed_kmer_store.h"

#define FILTER_ON 0

/* FIXME - Change this initial capacity using preprocessing step. */
static const int INITIAL_CAPACITY = 100000;

KmerStore::KmerStore(k_t k)
    : k(k),
    kmer_filter(INITIAL_CAPACITY, 0.001, 0,
            (bloom_filter_hash_func_t) kmer_hash_K),
    kmer_map(NULL)
{
    kmer_map = new HashMap<kmer_t, qual_counts_t>(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k));
}

void KmerStore::insert(qekmer_t* qekmer)
{
#if FILTER_ON
    if (!kmer_filter.check(qekmer->kmer)) {
        kmer_filter.add(qekmer->kmer);
    } else {
#endif
        kmer_map->try_insert(qekmer->kmer);
        if (qekmer->lqual > Q_MIN)
            kmer_map->map[qekmer->kmer].lquals[qekmer->exts.left]++;
        if (qekmer->rqual > Q_MIN)
            kmer_map->map[qekmer->kmer].rquals[qekmer->exts.right]++;
#if FILTER_ON
    }
#endif
}

void KmerStore::trim(TrimmedKmerStore& trimmed_kmer_store)
{
    for (map_type_t::iterator it = kmer_map->map.begin();
            it != kmer_map->map.end();
            it++) {
        kmer_t kmer = it->first;
        ext_map_t ext_map = it->second.ext_map(D_MIN);

        if (ext_map.valid())
            trimmed_kmer_store.insert(kmer, ext_map);
        else
            free(kmer);
    }
}

KmerStore::~KmerStore()
{
    kmer_map->map.clear();
    delete kmer_map;
}
