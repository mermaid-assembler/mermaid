#ifndef _CONTIG_H_
#define _CONTIG_H_

#include <cstdio>
#include <vector>

#include <boost/cstdint.hpp>

#include "kmer.h"

/* TODO - Vary this variable to make sure the subcontigs thing works */
const size_t SUBCONTIG_LEN = 4000;         /* Length in bases */

class Contig {
public:
    Contig(int32_t id);

    size_t get_len() { return len; }
    int32_t get_id() { return id; }
    int32_t get_next_id() { return next_id; }
    void set_next_id(int32_t next_id) { this->next_id = next_id; }

    void append_kmer(kmer_t kmer, k_t len);

    void fprint(FILE* outfile);
    void fprintln(FILE* outfile);

protected:
    std::vector<kmer_t> subcontigs;
    size_t len;
    int32_t id;
    int32_t next_id;
};

#endif /* _CONTIG_H_ */
