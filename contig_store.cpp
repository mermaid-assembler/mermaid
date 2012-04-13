#include <cstdio>
#include <vector>

#include "contig.h"
#include "contig_store.h"

using namespace std;

ContigStore::ContigStore()
    : contigs(), next_id(0)
{
}

Contig* ContigStore::get_new_contig()
{
    Contig* contig = new Contig(next_id++);
    contigs.push_back(contig);
    return contig;
}

void ContigStore::print_contigs(FILE* outfile)
{
    for (vector<Contig*>::iterator it = contigs.begin();
            it != contigs.end();
            it++) {
        fprintf(outfile, "id %d next_id %d len %lu\n", (*it)->get_id(), (*it)->get_next_id(), (*it)->get_len());
        (*it)->fprintln(outfile);
    }
}
