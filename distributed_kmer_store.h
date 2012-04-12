#ifndef _DISTRIBUTED_KMER_STORE_H_
#define _DISTRIBUTED_KMER_STORE_H_

#include <vector>

#include "kmer.h"
#include "hash_map.h"
#include "scalable_bloom_filter.h"
#include "contig.h"

typedef struct {
    count_t lquals[BASE::NUM_BASES];
    count_t rquals[BASE::NUM_BASES];
} qual_counts_t;

typedef struct {
    uint8_t left : 4;
    uint8_t right : 4;
} ext_map_t;

typedef struct {
    ext_map_t ext_map;
    int32_t contig_idx;
} kmer_info_t;

class DistributedKmerStore {
public:
    typedef HashMap<kmer_t, qual_counts_t>::map_type_t counts_map_type_t;;
    typedef HashMap<kmer_t, kmer_info_t>::map_type_t contig_map_type_t;;

    DistributedKmerStore(k_t k);

    /* Insert qekmer into the kmer hash map. */
    void insert(qekmer_t* qekmer);
    
    /* Trim kmers which don't appear at least D_MIN times. */
    void trim();

    void print_ufx(FILE* outfile);

    void build_contigs(std::vector<Contig*> contigs);

protected:
    void build_contig(int32_t contig_idx, Contig* contig, kmer_t beg_kmer, kmer_info_t& beg_kmer_info);

    k_t k;
    ScalableBloomFilter kmer_filter;
    HashMap<kmer_t, qual_counts_t>* counts_map;     /* Maps kmer to qual counts */
    HashMap<kmer_t, kmer_info_t>* contig_map;       /* Trimmed map that points
                                                       go contigs */
};

#endif /* _DISTRIBUTED_KMER_STORE_H_ */
