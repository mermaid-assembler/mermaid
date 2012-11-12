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
    ~ContigStore();

    void add(Contig* contig);
    //void add_to_final_contigs(Contig* contig);

    void fprint_contigs(FILE* outfile, size_t min_contig_len = Config::MIN_CONTIG_LEN);

    //typedef HashMap<kmer_t, Contig*>::map_type_t contig_map_type_t;
    //typedef contig_map_type_t::iterator iterator;

    ///* Find contig by its first k bases */
    //iterator find(kmer_t kmer, bool& revcmp_found);
    //iterator begin() { return contigs.map.begin(); };
    //iterator end() { return contigs.map.end(); };

    ///* This is necessary because we have two maps now... */
    //bool is_end(iterator it);
    //void trim();
    std::vector<Contig*> contigs;

protected:
    //// FIXME: Use a densehash instead of a sparsehash
    //HashMap<kmer_t, Contig*> contigs;
    ///* This hashes contigs based on the first k bases of their reverse
    // * complement.  Note that the contigs themselves are not duplicated, we
    // * just have two keys which map to each contig.*/
    //HashMap<kmer_t, Contig*> revmap_contigs;
    //std::vector<Contig*> final_contigs;

    k_t k;
};
#endif /* _CONTIG_STORE_H_ */
