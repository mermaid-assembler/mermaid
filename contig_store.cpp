#include <cassert>
#include <cstdio>
#include <vector>
#include <set>

#include "contig.h"
#include "contig_store.h"

using namespace std;

ContigStore::ContigStore(k_t k)
    : contigs(), k(k)
{
}

void ContigStore::add(Contig* contig)
{
    contigs.push_back(contig);
}

void ContigStore::fprint_contigs(FILE* outfile, size_t min_contig_len)
{
    for (vector<Contig*>::iterator it = contigs.begin();
            it != contigs.end();
            it++)
    {
        Contig* contig = *it;
        if (contig->s.size() >= min_contig_len)
            contig->fprint_fasta(outfile, Config::FASTA_TEXTWIDTH);
    }
}

ContigStore::~ContigStore()
{
    while (!contigs.empty()) {
        delete contigs.back();
        contigs.pop_back();
    }
}
