#include <cassert>
#include <cstdio>
#include <vector>
#include <set>

#include "contig.h"
#include "contig_store.h"

using namespace std;

static const int INITIAL_CAPACITY = 100000;
ContigStore::ContigStore(k_t k)
    : contigs(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K,
            kmer_size(k)), k(k)
{
}

void ContigStore::add_contig(Contig* contig)
{
    kmer_t kmer = (kmer_t) malloc(kmer_size(k));
    str2kmer(kmer, contig->s.c_str(), k);
    contigs.try_insert(kmer);
    contigs.map[kmer] = contig;
}

void ContigStore::print_contigs(FILE* outfile)
{
    for (contig_map_type_t::iterator it = contigs.map.begin();
            it != contigs.map.end();
            it++) {
        Contig* contig = it->second;

        if (contig->s.size() < MIN_CONTIG_LEN) continue;

        contig->fprint_fasta(outfile, FASTA_TEXTWIDTH);
    }
}

// vector versions:
#if 0
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
        if (contig->s.size() < MIN_CONTIG_LEN) continue;

        contig->fprint_fasta(outfile, FASTA_TEXTWIDTH);
    }
}
#endif
