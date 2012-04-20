#include <cassert>
#include <cstdio>
#include <vector>
#include <set>

#include "contig.h"
#include "contig_store.h"

using namespace std;

ContigStore::ContigStore()
    : contigs()
{
}

void ContigStore::add_contig(Contig* contig, int32_t next_id)
{
    /* TODO - Take out this assert for production code. */
    assert(contig->id == contigs.size());
    contigs.push_back(contig);

    if (next_id >= 0 && contig->id != (uint32_t) next_id) {
        Contig* next_contig = contigs[next_id];
        if (contig->extended_hash == next_contig->hash) {
           if (contig->can_join_contig(next_contig) &&
                   contig->join_contig(next_contig)) {
                delete next_contig;
                contigs[next_id] = NULL;
            }
        } else if (contig->extended_hash == next_contig->revcmp_hash) {
            next_contig->revcmp();
            if (contig->can_join_contig(next_contig) &&
                    contig->join_contig(next_contig)) {
                delete next_contig;
                contigs[next_id] = NULL;
            }
        }
    }
}

void ContigStore::print_contigs(FILE* outfile)
{
    for (vector<Contig*>::iterator it = contigs.begin();
            it != contigs.end();
            it++) {
        if (*it == NULL)
            continue;

        fprintf(outfile, "id %d len %lu\n", (*it)->id, (*it)->len);
        (*it)->fprintln(outfile);
    }
}
