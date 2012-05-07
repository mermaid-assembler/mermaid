#ifndef _KMER_EXT_MAP_H_
#define _KMER_EXT_MAP_H_

#include <cstdio>

#include "kmer.h"
#include "hash_map.h"

class ContigStore;

class KmerExtMap {
public:
    KmerExtMap(k_t k);
    ~KmerExtMap();

    void insert(kmer_t kmer, ext_map_t ext_map);
    void print_ufxs(FILE* outfile);
    void load_ufxs(FILE* infile);

    void build_contigs(ContigStore& contig_store);

protected:
    typedef HashMap<kmer_t, ext_map_t>::map_type_t map_type_t;

    void walk(Contig* contig);
    map_type_t::iterator lookup_kmer(kmer_t kmer, bool& used_revcmp);

    k_t k;
    HashMap<kmer_t, ext_map_t> hash_map;

    /* FIXME - Change this initial capacity using preprocessing step. */
    static const size_t INITIAL_CAPACITY = 100000;
};
#endif /* _KMER_EXT_MAP_H_ */
