#include "config.h"
#include "kmer.h"
#include "kmer_count_store.h"
#include "utils.h"

#define FILTER_ON 0

/* FIXME - Change this initial capacity using preprocessing step. */
static const int INITIAL_CAPACITY = 100000;

KmerCountStore::KmerCountStore(k_t k)
    : k(k),
    kmer_filter(INITIAL_CAPACITY, 0.001, 0,
            (bloom_filter_hash_func_t) kmer_hash_K),
    counts_map(NULL), contig_map(NULL)
{
    counts_map = new HashMap<kmer_t, qual_counts_t>(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k));
}

void KmerCountStore::insert(qekmer_t* qekmer)
{
#if FILTER_ON
    if (!kmer_filter.check(qekmer->kmer)) {
        kmer_filter.add(qekmer->kmer);
    } else {
#endif
        counts_map->try_insert(qekmer->kmer);
        if (qekmer->lqual > Q_MIN)
            counts_map->map[qekmer->kmer].lquals[qekmer->exts.left]++;
        if (qekmer->rqual > Q_MIN)
            counts_map->map[qekmer->kmer].rquals[qekmer->exts.right]++;
#if FILTER_ON
    }
#endif
}

/* Return a single byte which represents a bitmap for valid extensions. */
/* TODO - Explain this function more... */
static ext_map_t get_ext_map(qual_counts_t& qual_counts)
{
    ext_map_t ext_map = {0, 0};

    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        if (qual_counts.lquals[i] >= D_MIN)
            ext_map.left |= 1 << i;
        if (qual_counts.rquals[i] >= D_MIN)
            ext_map.right |= 1 << i;
    }
    return ext_map;
}

void KmerCountStore::trim()
{
    contig_map = new HashMap<kmer_t, kmer_info_t>(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k));

    for (counts_map_type_t::iterator it = counts_map->map.begin();
            it != counts_map->map.end();
            it++) {
        const kmer_t& kmer = it->first;
        qual_counts_t& qual_counts = it->second;

        ext_map_t ext_map = get_ext_map(qual_counts);
        if (ext_map.left || ext_map.right) {
            contig_map->map[kmer].ext_map = ext_map;
            contig_map->map[kmer].contig_idx = -1;
        } else {
            free(kmer);
        }
    }

    counts_map->map.clear();
    free(counts_map);
    counts_map = NULL;
}

void KmerCountStore::print_ufx(FILE* outfile)
{
    if (contig_map == NULL) {
        panic("You have not called trim on the kmer count store\n");
    }

    for (contig_map_type_t::iterator it = contig_map->map.begin();
            it != contig_map->map.end();
            it++) {
        const kmer_t& kmer = it->first;
        kmer_info_t& kmer_info = it->second;

        if (kmer_matches_str(kmer, "AAAAAACCTTACACACAGTGTTTTCTTTATTAGAAACTATT",
                    K)) {
            printf("%x %x\n", kmer_info.ext_map.left, kmer_info.ext_map.right);
        }
        char left_ext = 0;
        char right_ext = 0;

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (kmer_info.ext_map.left & 1 << i) {
                if (left_ext)
                    left_ext = 'F';
                else
                    left_ext = base2char((base) i);
            }
            if (kmer_info.ext_map.right & 1 << i) {
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
