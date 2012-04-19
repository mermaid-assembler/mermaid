#ifndef _CONTIG_STORE_H_
#define _CONTIG_STORE_H_

#include <cstdio>
#include <vector>

#include "contig.h"

class ContigStore {
public:
    ContigStore();

    void add_contig(Contig* contig);

    void print_contigs(FILE* outfile);

protected:
    std::vector<Contig*> contigs;
};
#endif /* _CONTIG_STORE_H_ */
