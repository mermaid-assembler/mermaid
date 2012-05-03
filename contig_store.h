#ifndef _CONTIG_STORE_H_
#define _CONTIG_STORE_H_

#include <cstdio>
#include <vector>
#include <map>

#include "contig.h"
#include "kmer.h"
#include "hash_map.h"

class ContigStore {
public:
    ContigStore(k_t k);

    void add_contig(Contig* contig);

    void print_contigs(FILE* outfile);

    size_t size() { return contigs.map.size(); };

    //typedef std::vector<Contig*>::iterator iterator;
    typedef HashMap<kmer_t, Contig*>::map_type_t contig_map_type_t;
    typedef contig_map_type_t::iterator iterator;

    /* Find contig by its first k bases */
    iterator find(kmer_t kmer, bool& revcmp_found);
    iterator begin() { return contigs.map.begin(); };
    iterator end() { return contigs.map.end(); };

protected:
    //std::vector<Contig*> contigs;
    HashMap<kmer_t, Contig*> contigs;

    k_t k;
};
#endif /* _CONTIG_STORE_H_ */
