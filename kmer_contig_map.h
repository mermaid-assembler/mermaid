#ifndef _KMER_CONTIG_MAP_H_
#define _KMER_CONTIG_MAP_H_

#include "config.h"
#include "kmer.h"
#include "contig.h"
#include "hash_map.h"
#include "contig_store.h"

class KmerContigMap {
public:
    KmerContigMap(k_t k);
    ~KmerContigMap();

    void insert(Contig* contig);
    void fprint_contigs(FILE* outfile, size_t min_contig_len = Config::MIN_CONTIG_LEN);

    void join_contigs(ContigStore& contig_store);

protected:
    typedef HashMap<kmer_t, Contig*>::map_type_t map_type_t;

    void walk(Contig* contig);
    map_type_t::iterator lookup_kmer(kmer_t kmer, bool& used_revcmp);

    k_t k;
    HashMap<kmer_t, Contig*> forward_map;
    HashMap<kmer_t, Contig*> revcmp_map;

    /* FIXME - Change this initial capacity using preprocessing step. */
    static const int INITIAL_CAPACITY = 100000;
};

#endif /* _KMER_CONTIG_MAP_H_ */
