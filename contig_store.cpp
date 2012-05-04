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
            kmer_size(k)),
      revmap_contigs(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K,
            kmer_size(k)),
      k(k)
{
}

void ContigStore::add_contig(Contig* contig)
{
    kmer_t kmer = (kmer_t) malloc(kmer_size(k));
    str2kmer(kmer, contig->s.c_str(), k);
    contigs.try_insert(kmer);
    contigs.map[kmer] = contig;

    kmer_t revcmp = (kmer_t) malloc(kmer_size(k));
    kmer_a tmp_kmer[kmer_size(k)];
    str2kmer(tmp_kmer, &contig->s.c_str()[contig->s.size()-k], k);
    revcmp_kmer(revcmp, tmp_kmer, k);
    revmap_contigs.try_insert(revcmp);
    revmap_contigs.map[revcmp] = contig;
}

void ContigStore::add_to_final_contigs(Contig* contig)
{
    final_contigs.push_back(contig);
}

void ContigStore::print_contigs(FILE* outfile, size_t min_contig_len)
{
    for (vector<Contig*>::iterator it = final_contigs.begin();
            it != final_contigs.end();
            it++)
    {
        Contig* contig = *it;
        if (contig->s.size() >= min_contig_len)
            contig->fprint_fasta(outfile, FASTA_TEXTWIDTH);
    }
}

ContigStore::iterator ContigStore::find(kmer_t kmer, bool& revcmp_found)
{
    iterator it;

    it = contigs.map.find(kmer);
    if (it != contigs.map.end()) {
        revcmp_found = false;
        return it;
    }

    it = revmap_contigs.map.find(kmer);
    if (it != revmap_contigs.map.end()) {
        revcmp_found = true;
        return it;
    }

    return it;
}

bool ContigStore::is_end(ContigStore::iterator it)
{
    return (it == contigs.map.end()) || (it == revmap_contigs.map.end());
}

void ContigStore::trim()
{
    for (contig_map_type_t::iterator it = contigs.map.begin();
            it != contigs.map.end();
            it++) {
        kmer_t kmer = it->first;
        Contig* contig = it->second;
        if (contig->s.size() == 0) delete contig;
        free(kmer);
    }
    for (contig_map_type_t::iterator it = revmap_contigs.map.begin();
            it != revmap_contigs.map.end();
            it++) {
        kmer_t kmer = it->first;
        free(kmer);
    }

    contigs.map.clear();
    revmap_contigs.map.clear();
}
