#include <cassert>
#include <cstdio>
#include <vector>

#include "contig.h"
#include "contig_store.h"

using namespace std;

ContigStore::ContigStore()
    : contigs()
{
}

void ContigStore::add_contig(Contig* contig)
{
    /* TODO - Take out this assert for production code. */
    assert(contig->id == (int32_t) contigs.size());
    contigs.push_back(contig);

    if (contig->next_id >= 0) {
        if (contig->id == contig->next_id) {
            contig->next_id = -1;
        } else {
            if (contig->join_contig(contigs[contig->next_id])) {
                delete contigs[contig->next_id];
                contigs[contig->next_id] = NULL;
                contig->next_id = -1;
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

        fprintf(outfile, "id %d next_id %d len %lu\n", (*it)->id, (*it)->next_id, (*it)->len);
        (*it)->fprintln(outfile);
    }
}
