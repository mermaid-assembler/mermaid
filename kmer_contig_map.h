#ifndef _KMER_CONTIG_MAP_H_
#define _KMER_CONTIG_MAP_H_

#include "config.h"
#include "kmer.h"
#include "contig.h"
#include "hash_map.h"

class KmerContigMap {
public:
    KmerContigMap(k_t k);

    void insert(Contig* contig);
    void fprint_contigs(FILE* outfile, size_t min_contig_len = MIN_CONTIG_LEN);

protected:
    typedef HashMap<kmer_t, Contig*>::map_type_t map_type_t;

    k_t k;
    HashMap<kmer_t, Contig*>* hash_map;

    /* FIXME - Change this initial capacity using preprocessing step. */
    static const int INITIAL_CAPACITY = 100000;
};

#endif /* _KMER_CONTIG_MAP_H_ */
