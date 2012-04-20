#ifndef _CONTIG_H_
#define _CONTIG_H_

#include <cstdio>
#include <vector>

#include <boost/cstdint.hpp>

#include "kmer.h"

/* WARNING - This must be larger than k. */
const size_t SUBCONTIG_LEN = 400;         /* Length in bases */

class Contig {
public:
    static void set_k(k_t k) { Contig::k = k; }

    Contig();
    ~Contig();

    void append_base(base b);
    void append_first_kmer(kmer_t kmer);
    /* Appends a contig to this contig. Returns true if the join succeeded;
     * false otherwise. The join may fail if next_contig's left extension
     * doesn't match up.
     */
    bool can_join_contig(Contig* next_contig);
    bool join_contig(Contig* next_contig);

    /* Reverse complements the current contig. */
    void revcmp(void);

    void generate_hashes(void);

    void fprint(FILE* outfile);
    void fprintln(FILE* outfile);

public:
    exts_t exts;
    size_t len;
    uint32_t id;
    size_t hash;
    size_t extended_hash;       /* Hash of kmer when extended to right. */
    size_t revcmp_hash;

protected:
    void get_extended_kmer(kmer_t extended_kmer);

    static uint32_t id_generator;
    static k_t k;
    static size_t seed;

    std::vector<kmer_t> subcontigs;
};

#endif /* _CONTIG_H_ */
