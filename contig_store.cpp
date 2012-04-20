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
    next_id = -1;
    /* TODO - Take out this assert for production code. */
    assert(contig->id == contigs.size());
    contigs.push_back(contig);

    if (next_id >= 0 && contig->id != (uint32_t) next_id) {
        Contig* next_contig = contigs[next_id];
        if (contig->extended_hash == next_contig->hash) {
           if (contig->can_join_contig(next_contig) &&
                   contig->join_contig(next_contig)) {
                delete next_contig;
                contigs[next_id] = contig;
            }
        } else if (contig->extended_hash == next_contig->revcmp_hash) {
            next_contig->revcmp();
            if (contig->can_join_contig(next_contig) &&
                    contig->join_contig(next_contig)) {
                delete next_contig;
                contigs[next_id] = contig;
            }
        }
    }
}

void ContigStore::print_contigs(FILE* outfile)
{
    set<uint32_t> printed_ids;

    for (vector<Contig*>::iterator it = contigs.begin();
            it != contigs.end();
            it++) {
        if (printed_ids.find((*it)->id) != printed_ids.end())
            continue;

        fprintf(outfile, "id %d len %lu\n", (*it)->id, (*it)->len);
        (*it)->fprintln(outfile);
        printed_ids.insert((*it)->id);
    }
}
