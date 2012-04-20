#ifndef _CONTIG_STORE_H_
#define _CONTIG_STORE_H_

#include <cstdio>
#include <vector>
#include <map>

#include "contig.h"

class ContigStore {
public:
    ContigStore();

    void add_contig(Contig* contig, int32_t next_id);
    void join_contigs(void);

    void print_contigs(FILE* outfile);

protected:
    std::map<uint32_t, Contig*> contigs;
    std::map<size_t, std::map<uint32_t, Contig*> > hash_id_map;
};
#endif /* _CONTIG_STORE_H_ */
