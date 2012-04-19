#ifndef _CONTIG_H_
#define _CONTIG_H_

#include <cstdio>
#include <vector>

#include <boost/cstdint.hpp>

#include "kmer.h"

const size_t SUBCONTIG_LEN = 400;         /* Length in bases */

class Contig {
public:
    static void set_k(k_t k) { Contig::k = k; }

    Contig();

    void append_base(base b);
    void append_first_kmer(kmer_t kmer);
    /* Appends a contig to this contig. Returns true if the join succeeded;
     * false otherwise. The join may fail if next_contig's left extension
     * doesn't match up.
     */
    bool join_contig(Contig* next_contig);

    void fprint(FILE* outfile);
    void fprintln(FILE* outfile);

public:
    exts_t exts;
    size_t len;
    int32_t id;
    int32_t next_id;

protected:
    static int32_t id_generator;
    static k_t k;

    std::vector<kmer_t> subcontigs;
};

#endif /* _CONTIG_H_ */
