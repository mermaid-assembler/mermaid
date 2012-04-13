#include "config.h"
#include "kmer.h"
#include "kmer_count_store.h"
#include "utils.h"
#include "contig_store.h"

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
            contig_map->map[kmer].contig_id = -1;
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

static bool can_use_in_contig(kmer_info_t& kmer_info)
{
    uint8_t valid_left_bases = 0;
    uint8_t valid_right_bases = 0;

    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        if (kmer_info.ext_map.left & 1 << i)
            valid_left_bases++;
        if (kmer_info.ext_map.right & 1 << i)
            valid_right_bases++;
    }

    return valid_left_bases == 1 && valid_right_bases == 1;
}

/* WARNING - This function assumes the given side of the ext_map has only one
 * bit set. */
static base ext_map_side2base(uint8_t side)
{
    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        if (side & 1 << i)
            return (base) i;
    }

    panic("Given ext_map side didn't have any bits set in %s\n", __func__);
}

void KmerCountStore::build_contig(Contig* contig, kmer_t beg_kmer, kmer_info_t& beg_kmer_info)
{
    if (!can_use_in_contig(beg_kmer_info))
        return;

    kmer_a subcontig[SUBCONTIG_LEN];
    kmer_a tmp_kmer[kmer_size(k)];
    size_t idx;
    base b;

    set_base(subcontig, 0, (base) ext_map_side2base(beg_kmer_info.ext_map.left));
    for_base_in_kmer(b, beg_kmer, k) {
        set_base(subcontig, 1 + b_i_, b);
    } end_for;
    contig->append_kmer(subcontig, k + 1);
    beg_kmer_info.contig_id = contig->get_id();
    idx = k + 1;

    kmer_info_t& kmer_info = beg_kmer_info;
    while (1) {
        set_base(subcontig, idx, ext_map_side2base(kmer_info.ext_map.right));
        idx++;
        for_base_in_kmer_from(b, subcontig, k, idx - k) {
            set_base(tmp_kmer, b_i_, b);
        } end_for;
        contig->append_kmer(tmp_kmer, k);

        contig_map_type_t::iterator it = contig_map->map.find(tmp_kmer);
        if (it == contig_map->map.end())
            return;

        kmer_info = it->second;

        if (!can_use_in_contig(kmer_info))
            return;
        /* If not reflexive */
        if (ext_map_side2base(kmer_info.ext_map.left) !=
                get_base(subcontig, idx - (k + 1)))
            return;

        if (kmer_info.contig_id >= 0) {
            contig->set_next_id(kmer_info.contig_id);
            return;
        } else {
            kmer_info.contig_id = contig->get_id();
        }

        if (idx == SUBCONTIG_LEN) {
            memcpy(subcontig, tmp_kmer, kmer_size(k));
            idx = 0;
        }
    }
}

void KmerCountStore::build_contigs(ContigStore& contig_store)
{
    for (contig_map_type_t::iterator it = contig_map->map.begin();
            it != contig_map->map.end();
            it++) {

        kmer_t kmer = it->first;
        kmer_info_t& kmer_info = it->second;

        if (!can_use_in_contig(kmer_info))
            continue;
        if (kmer_info.contig_id >= 0)
            continue;

        Contig* contig = contig_store.get_new_contig();
        build_contig(contig, kmer, kmer_info);
    }
}

/* TODO - Version of build_contig which doesn't constantly copy the next kmer
 * into a temporary buffer. */
#if 0
void KmerCountStore::build_contig(int32_t contig_idx, Contig* contig, const kmer_t& beg_kmer, kmer_info_t& beg_kmer_info)
{
    if (!can_use_in_contig(beg_kmer_info))
        return;

    /* We don't want to keep copying the next kmer to a temporary buffer and
     * using that to index into the map. Instead, we want an array which will
     * already have the previous (k-1) bases filled out. We want BASES_PER_BYTE
     * of these arrays because kmers start at different byte offsets depending
     * on the parity.
     * FIXME - Explain better...
     */
    kmer_a subcontigs[BASES_PER_BYTE][SUBCONTIG_LEN];
    size_t idx;
    uint8_t parity = 0;
    base b;

    set_base(subcontigs[0], 0, (base) ext_map_side2base(beg_kmer_info.ext_map.left));
    for_base_in_kmer(b, beg_kmer, k) {
        set_base(subcontigs[0], 1 + b_i_, b);
    } end_for;
    idx = k + 1;
    contig->append_kmer(beg_kmer, k + 1);
    beg_kmer_info.contig_idx = contig_idx;
    parity = (parity + 1) % BASES_PER_BYTE;

    for (uint8_t i = 1; i < BASES_PER_BYTE; i++) {
        for_base_in_kmer_from(b, beg_kmer, k - (i - 1), i - 1) {
            set_base(subcontigs[i], b_i_, b);
        } end_for;
    }

    const kmer_t& kmer = beg_kmer;
    kmer_info_t& kmer_info = beg_kmer_info;
    while (1) {
        for (uint8_t i = 0; i < BASES_PER_BYTE; i++) {
            b = (base) ext_map_side2base(kmer_info.ext_map.right);
            set_base(subcontigs[i], idx - i, b);
        }
        kmer_t next_kmer = &subcontigs[parity][idx

        contig_map_type_t::iterator it = contig_map.find(parity
    }
}
#endif