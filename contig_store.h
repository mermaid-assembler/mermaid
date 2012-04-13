#ifndef _CONTIG_STORE_H_
#define _CONTIG_STORE_H_

#include <vector>

#include "contig.h"

class ContigStore {
public:
    ContigStore();

    Contig* get_new_contig();

protected:
    std::vector<Contig*> contigs;
    int32_t next_id;
};
#endif /* _CONTIG_STORE_H_ */
