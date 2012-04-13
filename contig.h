#ifndef _CONTIG_H_
#define _CONTIG_H_

#include <cstdio>
#include <vector>

#include <boost/cstdint.hpp>

#include "kmer.h"

const size_t SUBCONTIG_LEN = 400;         /* Length in bases */

class Contig {
public:
    Contig(int32_t id);

    void append_kmer(kmer_t kmer, k_t len);

    void fprint(FILE* outfile);
    void fprintln(FILE* outfile);

protected:
    std::vector<kmer_t> subcontigs;

public:
    struct {
        uint8_t left : 4;
        uint8_t right : 4;
    } exts;
    size_t len;
    int32_t id;
    int32_t next_id;
};

#endif /* _CONTIG_H_ */
