#include <cassert>
#include <cstdio>
#include <vector>
#include <set>

#include "contig.h"
#include "contig_store.h"

using namespace std;

ContigStore::ContigStore()
    : contigs(), hash_id_map()
{
}

void ContigStore::add_contig(Contig* contig, int32_t next_id)
{
    /* TODO - Take out this assert for production code. */
    contigs[contig->id] = contig;

    if (next_id >= 0 && contig->id != (uint32_t) next_id) {
        Contig* next_contig = contigs[next_id];
        if (contig->extended_hash == next_contig->hash) {
           if (contig->can_join_contig(next_contig) &&
                   contig->join_contig(next_contig)) {
                delete next_contig;
                contigs.erase(next_id);
            }
        } else if (contig->extended_hash == next_contig->revcmp_hash) {
            next_contig->revcmp();
            if (contig->can_join_contig(next_contig) &&
                    contig->join_contig(next_contig)) {
                delete next_contig;
                contigs.erase(next_id);
            }
        }
    }

    hash_id_map[contig->hash][contig->id] = contig;
    hash_id_map[contig->revcmp_hash][contig->id] = contig;
}

void ContigStore::join_contigs(void)
{
}

void ContigStore::print_contigs(FILE* outfile)
{
    for (map<uint32_t, Contig*>::iterator it = contigs.begin();
            it != contigs.end();
            it++) {
        Contig* contig = it->second;

        fprintf(outfile, "id %d len %lu\n", contig->id, contig->len);
        contig->fprintln(outfile);
    }
}
