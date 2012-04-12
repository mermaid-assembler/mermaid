#ifndef _CONTIG_H_
#define _CONTIG_H_

#include <vector>

#include <boost/cstdint.hpp>

#include "kmer.h"

const size_t SUBCONTIG_LEN = 4000;         /* Length in bases */

class Contig {
public:
    Contig();

    void append_kmer(kmer_t kmer, k_t len);

protected:
    std::vector<kmer_t> subcontigs;
    size_t len;
    int32_t next_contig_id;
};

#endif /* _CONTIG_H_ */
