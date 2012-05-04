#include "kmer_contig_map.h"

KmerContigMap::KmerContigMap(k_t k)
    : k(k), forward_map(NULL), revcmp_map(NULL)
{
    forward_map = new HashMap<kmer_t, Contig*>(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k));
    revcmp_map = new HashMap<kmer_t, Contig*>(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k));
}

void KmerContigMap::insert(Contig* contig)
{
    kmer_t kmer = (kmer_t) malloc(kmer_size(k));
    for (size_t i = 0; i < k; i++)
        set_base(kmer, i, char2base(contig->s[i]));
    forward_map->map[kmer] = contig;
    contig->revcmp();
    kmer_t revcmp = (kmer_t) malloc(kmer_size(k));
    for (size_t i = 0; i < k; i++)
        set_base(revcmp, i, char2base(contig->s[i]));
    revcmp_map->map[revcmp] = contig;
}

void KmerContigMap::fprint_contigs(FILE* outfile, size_t min_contig_len)
{
    for (map_type_t::iterator it = forward_map->map.begin();
            it != forward_map->map.end();
            it++)
    {
        Contig* contig = it->second;
        if (contig->s.size() >= min_contig_len)
            contig->fprint_fasta(outfile, FASTA_TEXTWIDTH);
    }
}

void KmerContigMap::join_contigs(ContigStore& contig_store)
{
    for (map_type_t::iterator it = forward_map->map.begin();
            it != forward_map->map.end();
            it++) {
        Contig* contig = it->second;

        if (contig->s.size() == 0) continue;

        walk(contig);
        contig->revcmp();
        walk(contig);
        contig_store.add(contig);
    }
}

KmerContigMap::map_type_t::iterator KmerContigMap::lookup_kmer(kmer_t kmer, bool& used_revcmp)
{
    map_type_t::iterator it;

    it = forward_map->map.find(kmer);
    if (it != forward_map->map.end()) {
        used_revcmp = false;
        return it;
    }

    it = revcmp_map->map.find(kmer);
    if (it != revcmp_map->map.end()) {
        used_revcmp = true;
        return it;
    }

    return it;
}

void KmerContigMap::walk(Contig* contig)
{
    kmer_a kmer[kmer_size(k)];
    bool used_revcmp;
    base left_ext;

    while (1) {
        contig->get_ext_kmer(kmer);
        map_type_t::iterator it = lookup_kmer(kmer, used_revcmp);
        if (!used_revcmp) {
            if (it == forward_map->map.end()) break;
        } else {
            if (it == revcmp_map->map.end()) break;
        }

        Contig*& next_contig = it->second;
        if (next_contig == contig || next_contig->s.size() == 0) break;

        if (!used_revcmp) {
            left_ext = next_contig->left_ext;
        } else {
            left_ext = inv_base(next_contig->right_ext);
        }

        if (!contig->check_next_left_ext(left_ext)) break;

        if (used_revcmp)
            next_contig->revcmp();

        contig->append(next_contig);
        next_contig->s.clear();
    }
}
