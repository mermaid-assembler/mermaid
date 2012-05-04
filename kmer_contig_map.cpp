#include "kmer_contig_map.h"

KmerContigMap::KmerContigMap(k_t k)
    : k(k), hash_map(NULL)
{
    hash_map = new HashMap<kmer_t, Contig*>(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k));
}

void KmerContigMap::insert(Contig* contig)
{
    kmer_t kmer = (kmer_t) malloc(kmer_size(k));

    for (size_t i = 0; i < k; i++)
        set_base(kmer, i, char2base(contig->s[i]));

    hash_map->map[kmer] = contig;
}

void KmerContigMap::fprint_contigs(FILE* outfile, size_t min_contig_len)
{
    for (map_type_t::iterator it = hash_map->map.begin();
            it != hash_map->map.end();
            it++)
    {
        Contig* contig = it->second;
        if (contig->s.size() >= min_contig_len)
            contig->fprint_fasta(outfile, FASTA_TEXTWIDTH);
    }
}
