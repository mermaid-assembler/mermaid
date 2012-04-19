#include <cstdio>
#include <vector>

#include "contig.h"
#include "contig_store.h"

using namespace std;

ContigStore::ContigStore()
    : contigs()
{
}

Contig* ContigStore::get_new_contig()
{
    Contig* contig = new Contig();
    contigs.push_back(contig);
    return contig;
}

void ContigStore::print_contigs(FILE* outfile)
{
    for (vector<Contig*>::iterator it = contigs.begin();
            it != contigs.end();
            it++) {
        fprintf(outfile, "id %d next_id %d len %lu\n", (*it)->id, (*it)->next_id, (*it)->len);
        (*it)->fprintln(outfile);
    }
}
