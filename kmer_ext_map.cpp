#include "config.h"
#include "contig.h"
#include "kmer_ext_map.h"
#include "contig_store.h"

KmerExtMap::KmerExtMap(k_t k)
    :k(k), hash_map(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k))
{
}

void KmerExtMap::insert(kmer_t kmer, ext_map_t ext_map)
{
    hash_map.map[kmer] = ext_map;
}

void KmerExtMap::print_ufxs(FILE* outfile)
{
    for (map_type_t::iterator it = hash_map.map.begin();
            it != hash_map.map.end();
            it++) {
        kmer_t kmer = it->first;
        ext_map_t ext_map = it->second;

        char left_ext = 0;
        char right_ext = 0;

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (ext_map.left & 1 << i) {
                if (left_ext)
                    left_ext = 'F';
                else
                    left_ext = base2char((base) i);
            }
            if (ext_map.right & 1 << i) {
                if (right_ext)
                    right_ext = 'F';
                else
                    right_ext = base2char((base) i);
            }
        }
        if (!left_ext)
            left_ext = 'X';
        if (!right_ext)
            right_ext = 'X';

        fprint_kmer(outfile, kmer, k);
        fprintf(outfile, "\t%c%c\n", left_ext, right_ext);
        kmer_a revcmp[kmer_size(k)];
        revcmp_kmer(revcmp, kmer, k);
        fprint_kmer(outfile, revcmp, k);
        fprintf(outfile, "\t");
        if (right_ext == 'F' || right_ext == 'X')
            fprintf(outfile, "%c", right_ext);
        else
            fprintf(outfile, "%c", inv_base(right_ext));
        if (left_ext == 'F' || left_ext == 'X')
            fprintf(outfile, "%c", left_ext);
        else
            fprintf(outfile, "%c", inv_base(left_ext));
        fprintf(outfile, "\n");
    }
}

void KmerExtMap::load_ufxs(FILE* infile)
{
    char kmer_str[k + 1];
    char left_ext;
    char right_ext;
    kmer_a kmer[kmer_size(k)];

    while (fscanf(infile, "%s %c%c", kmer_str, &left_ext, &right_ext) != EOF) {
        str2kmer(kmer, kmer_str, k);
        if (!is_canonical_kmer(kmer, k)) continue;

        ext_map_t ext_map;
        ext_map.left = valid_base(left_ext) ? 1 << (uint8_t) char2base(left_ext) : 0;
        ext_map.right = valid_base(right_ext) ? 1 << (uint8_t) char2base(right_ext) : 0;
        hash_map.try_insert(kmer);
        hash_map.map[kmer] = ext_map;
    }
}

void KmerExtMap::build_contigs(ContigStore& contig_store)
{
    for (map_type_t::iterator it = hash_map.map.begin();
            it != hash_map.map.end();
            it++) {
        kmer_t kmer = it->first;
        ext_map_t ext_map = it->second;

        if (!ext_map.valid() || !ext_map.is_uu()) continue;

        Contig* contig = new Contig(kmer);
        contig->left_ext = ext_map.left_ext();
        contig->right_ext = ext_map.right_ext();
        walk(contig);
        contig->revcmp();
        walk(contig);
        contig_store.add(contig);
    }
}

KmerExtMap::map_type_t::iterator KmerExtMap::lookup_kmer(kmer_t kmer, bool& used_revcmp)
{
    map_type_t::iterator it;
    kmer_a revcmp[kmer_size(k)];

    it = hash_map.map.find(kmer);
    if (it != hash_map.map.end()) {
        used_revcmp = false;
        return it;
    }

    revcmp_kmer(revcmp, kmer, k);
    it = hash_map.map.find(revcmp);
    if (it != hash_map.map.end()) {
        used_revcmp = true;
        return it;
    }

    return it;
}

void KmerExtMap::walk(Contig* contig)
{
    kmer_a kmer[kmer_size(k)];
    bool used_revcmp;
    base left_ext;
    base right_ext;

    while (1) {
        contig->get_ext_kmer(kmer);
        map_type_t::iterator it = lookup_kmer(kmer, used_revcmp);
        if (it == hash_map.map.end()) break;

        ext_map_t& ext_map = it->second;
        if (!ext_map.valid() || !ext_map.is_uu()) break;

        if (!used_revcmp) {
            left_ext = ext_map.left_ext();
            right_ext = ext_map.right_ext();
        } else {
            left_ext = inv_base(ext_map.right_ext());
            right_ext = inv_base(ext_map.left_ext());
        }

        if (!contig->check_next_left_ext(left_ext)) break;

        contig->s += base2char(contig->right_ext);
        ext_map.invalidate();
        contig->right_ext = right_ext;
    }
}


KmerExtMap::~KmerExtMap()
{
    for (map_type_t::iterator it = hash_map.map.begin();
            it != hash_map.map.end();
            it++) {
        kmer_t kmer = it->first;
        free(kmer);
    }
}
