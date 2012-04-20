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

void ContigStore::add_contig(Contig* contig)
{
    contigs.push_back(contig);
}

void ContigStore::print_contigs(FILE* outfile)
{
    for (vector<Contig*>::iterator it = contigs.begin();
            it != contigs.end();
            it++) {
        Contig* contig = *it;

        fprintf(outfile, "id %d len %lu\n", contig->id, contig->s.size());
        contig->fprintln(outfile);
    }
}
