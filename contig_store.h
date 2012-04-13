#ifndef _CONTIG_STORE_H_
#define _CONTIG_STORE_H_

#include <cstdio>
#include <vector>

#include "contig.h"

class ContigStore {
public:
    ContigStore();

    Contig* get_new_contig();

    void print_contigs(FILE* outfile);

protected:
    std::vector<Contig*> contigs;
    int32_t next_id;
};
#endif /* _CONTIG_STORE_H_ */
