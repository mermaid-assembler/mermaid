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
