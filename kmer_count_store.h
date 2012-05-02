#ifndef _KMER_COUNT_STORE_H_
#define _KMER_COUNT_STORE_H_

#include <vector>

#include "kmer.h"
#include "hash_map.h"
#include "scalable_bloom_filter.h"
#include "contig.h"
#include "contig_store.h"

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
    bool contig_found;
} kmer_info_t;

typedef struct {
    int32_t next_id;
    bool next_revcmp;
} next_contig_info_t;

class KmerCountStore {
public:
    typedef HashMap<kmer_t, qual_counts_t>::map_type_t counts_map_type_t;
    typedef HashMap<kmer_t, kmer_info_t>::map_type_t contig_map_type_t;

    KmerCountStore(k_t k);

    /* Insert qekmer into the kmer hash map. */
    void insert(qekmer_t* qekmer);
    
    /* Trim kmers which don't appear at least D_MIN times. */
    void trim();

    void print_ufxs(FILE* outfile);
    void load_ufxs(FILE* infile);

    void build_contigs(ContigStore& contig_store);

protected:
    /* Returns id of the next contig. */
    void walk(Contig* contig, base next_ext);
    /* Takes the kmer from scratch_kmer and looks for the next kmer in
     * contig_map. If found, returns true; otherwise false is returned.
     * revcmp_found represents whether the original kmer was found or whether
     * reverse complement was found. kmer_info holds the information from the
     * next kmer. If the next kmer could not be found, revcmp_found and
     * kmer_info should not be disturbed. */
    contig_map_type_t::iterator get_next_kmer(kmer_t kmer, bool& revcmp_found);

    k_t k;
    ScalableBloomFilter kmer_filter;
    HashMap<kmer_t, qual_counts_t>* counts_map;     /* Maps kmer to qual counts */
    HashMap<kmer_t, kmer_info_t>* contig_map;       /* Trimmed map that points
                                                       go contigs */
};

#endif /* _KMER_COUNT_STORE_H_ */
