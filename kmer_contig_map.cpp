#include "kmer_contig_map.h"

KmerContigMap::KmerContigMap(k_t k)
    : k(k),
    forward_map(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k)),
    revcmp_map(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k))
{
}

void KmerContigMap::insert(Contig* contig)
{
    kmer_t kmer = (kmer_t) malloc(kmer_size(k));
    str2kmer(kmer, contig->s.c_str(), k);
    forward_map.map[kmer] = contig;

    kmer_t revcmp = (kmer_t) malloc(kmer_size(k));
    kmer_a tmp_kmer[kmer_size(k)];
    str2kmer(tmp_kmer, &contig->s.c_str()[contig->s.size()-k], k);
    revcmp_kmer(revcmp, tmp_kmer, k);
    revcmp_map.map[revcmp] = contig;
}

void KmerContigMap::fprint_contigs(FILE* outfile, size_t min_contig_len)
{
    for (map_type_t::iterator it = forward_map.map.begin();
            it != forward_map.map.end();
            it++)
    {
        Contig* contig = it->second;
        if (contig->s.size() >= min_contig_len)
            contig->fprint_fasta(outfile, Config::FASTA_TEXTWIDTH);
    }
}

void KmerContigMap::join_contigs(ContigStore& contig_store)
{
    for (map_type_t::iterator it = forward_map.map.begin();
            it != forward_map.map.end();
            it++) {
        if (it->second->s.size() == 0) continue;

        Contig* contig = new Contig(it->second);
        walk(contig);
        contig->revcmp();
        walk(contig);
        contig_store.add(contig);

        it->second->s.clear();
    }
}

KmerContigMap::map_type_t::iterator KmerContigMap::lookup_kmer(kmer_t kmer, bool& used_revcmp)
{
    map_type_t::iterator it;

    it = forward_map.map.find(kmer);
    if (it != forward_map.map.end()) {
        used_revcmp = false;
        return it;
    }

    it = revcmp_map.map.find(kmer);
    if (it != revcmp_map.map.end()) {
        used_revcmp = true;
        return it;
    }

    return it;
}

void KmerContigMap::walk(Contig* contig)
{
    kmer_a kmer[kmer_size(k)];
    bool used_revcmp;

    while (1) {
        contig->get_ext_kmer(kmer);
        map_type_t::iterator it = lookup_kmer(kmer, used_revcmp);
        if (it == forward_map.map.end() || it == revcmp_map.map.end())
            break;

        Contig*& next_contig = it->second;
        if (next_contig == contig || next_contig->s.size() == 0) break;

        if (used_revcmp)
            next_contig->revcmp();

        if (!contig->check_next_left_ext(next_contig->left_ext)) {
            if (used_revcmp)
                next_contig->revcmp();
            break;
        }

        contig->append(next_contig);
        next_contig->s.clear();
    }
}

KmerContigMap::~KmerContigMap()
{
    for (map_type_t::iterator it = forward_map.map.begin();
            it != forward_map.map.end();
            it++) {
        kmer_t kmer = it->first;
        Contig* contig = it->second;
        free(kmer);
        if (contig->s.size() == 0)
            delete contig;
    }

    for (map_type_t::iterator it = revcmp_map.map.begin();
            it != revcmp_map.map.end();
            it++) {
        kmer_t kmer = it->first;
        free(kmer);
    }
}
