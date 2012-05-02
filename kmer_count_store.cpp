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
            contig_map->map[kmer].contig_found = false;
        } else {
            free(kmer);
        }
    }

    counts_map->map.clear();
    delete counts_map;
    counts_map = NULL;
}

void KmerCountStore::print_ufxs(FILE* outfile)
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

static ext_map_t base2ext_map(char l, char r)
{
    ext_map_t ext_map = {0, 0};
    switch (l)
    {
        case 'A': ext_map.left = 1; break;
        case 'C': ext_map.left = 2; break;
        case 'G': ext_map.left = 4; break;
        case 'T': ext_map.left = 8; break;
        case 'F': ext_map.left = 15; break; // Not sure what I should do with this
        case 'X': ext_map.left = 0; break;
        default: panic ("Invalid base");
    }
    switch (r)
    {
        case 'A': ext_map.right = 1; break;
        case 'C': ext_map.right = 2; break;
        case 'G': ext_map.right = 4; break;
        case 'T': ext_map.right = 8; break;
        case 'F': ext_map.right = 15; break; // Not sure what I should do with this
        case 'X': ext_map.right = 0; break;
        default: panic ("Invalid base");
    }
    return ext_map;
}

bool is_canonical_kmer(kmer_t kmer, k_t k)
{
    kmer_a revcmp_a[k];
    kmer_t revcmp = revcmp_a;
    revcmp_kmer(revcmp, kmer, k);
    return cmp_kmer(kmer, revcmp,k) <= 0;
}

/* Move this some place neater */
void KmerCountStore::load_ufxs(FILE* infile)
{
    char kmer_str[k + 1];
    char left_ext;
    char right_ext;

    contig_map = new HashMap<kmer_t, kmer_info_t>(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K,
            (hash_map_eq_func_t) kmer_eq_K, kmer_size(k));
    while (fscanf(infile, "%s %c%c", kmer_str, &left_ext, &right_ext) != EOF) {
        kmer_t kmer = (kmer_t) malloc(kmer_size(k));
        str2kmer(kmer, kmer_str, k);
        if (!is_canonical_kmer(kmer, k)) continue;
        ext_map_t ext_map = base2ext_map(left_ext, right_ext);
        contig_map->map[kmer].ext_map = ext_map;
        contig_map->map[kmer].contig_found  = false;
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
static base ext_map_side_to_base(uint8_t side)
{
    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        if (side & 1 << i)
            return (base) i;
    }

    panic("Given ext_map side didn't have any bits set in %s\n", __func__);
}

KmerCountStore::contig_map_type_t::iterator KmerCountStore::get_next_kmer(kmer_t kmer, bool& revcmp_found)
{
    contig_map_type_t::iterator it;
    kmer_a revcmp[kmer_size(k)];

    it = contig_map->map.find(kmer);
    if (it != contig_map->map.end()) {
        revcmp_found = false;
        return it;
    }

    revcmp_kmer(revcmp, kmer, k);
    it = contig_map->map.find(revcmp);
    if (it != contig_map->map.end()) {
        revcmp_found = true;
        return it;
    }

    return it;
}

void KmerCountStore::walk(Contig* contig, base next_ext)
{
    kmer_a kmer[kmer_size(k)];
    bool revcmp_found;
    base left_ext;
    base right_ext;

    while (1) {
        for (size_t i = 0; i < k - 1; i++) {
            base b = char2base(contig->s[contig->s.size() - (k - 1) + i]);
            set_base(kmer, i, b);
        }
        set_base(kmer, k - 1, next_ext);

        contig_map_type_t::iterator it = get_next_kmer(kmer, revcmp_found);
        if (it == contig_map->map.end()) break;

        kmer_info_t& kmer_info = it->second;
        if (!can_use_in_contig(kmer_info)) break;
        if (kmer_info.contig_found) break;

        if (!revcmp_found) {
            left_ext = ext_map_side_to_base(kmer_info.ext_map.left);
            right_ext = ext_map_side_to_base(kmer_info.ext_map.right);
        } else {
            left_ext = inv_base(ext_map_side_to_base(kmer_info.ext_map.right));
            right_ext = inv_base(ext_map_side_to_base(kmer_info.ext_map.left));
        }

        if (!contig->check_next_left_ext(left_ext)) break;

        contig->s += base2char(next_ext);
        kmer_info.contig_found = true;
        next_ext = right_ext;
    }
    contig->right_ext = next_ext;
}

#if 0
int32_t KmerCountStore::build_contig(Contig* contig, kmer_t beg_kmer, kmer_info_t& beg_kmer_info)
{
    kmer_a subcontig[kmer_size(SUBCONTIG_LEN)];
    kmer_a cur_kmer[kmer_size(k)];
    exts_t exts;
    size_t idx;
    base b;
    bool revcmp_found;
    int32_t next_id = -1;

    contig->exts.left = ext_map_side2base(beg_kmer_info.ext_map.left);
    exts.left = contig->exts.left;
    for_base_in_kmer(b, beg_kmer, k) {
        set_base(subcontig, b_i_, b);
    } end_for;
    contig->append_first_kmer(subcontig);
    beg_kmer_info.contig_id = contig->id;
    idx = k;
    exts.right = ext_map_side2base(beg_kmer_info.ext_map.right);

    while (1) {
        set_base(subcontig, idx, exts.right);
        idx++;
        for_base_in_kmer_from(b, subcontig, k, idx - k) {
            set_base(cur_kmer, b_i_, b);
        } end_for;

        contig_map_type_t::iterator it = get_next_kmer(cur_kmer, revcmp_found);
        if (it == contig_map->map.end())
            break;
           
        kmer_info_t& cur_kmer_info = it->second;
        if (!can_use_in_contig(cur_kmer_info))
            break;

        if (!revcmp_found)
            exts.left = ext_map_side2base(cur_kmer_info.ext_map.left);
        else
            exts.left = inv_base(ext_map_side2base(cur_kmer_info.ext_map.right));
        
        if (get_base(subcontig, idx - k - 1) != exts.left)
            break;

        if (cur_kmer_info.contig_id >= 0) {
            next_id = cur_kmer_info.contig_id;
            break;
        } else {
            cur_kmer_info.contig_id = contig->id;
        }

        if (!revcmp_found)
            exts.right = ext_map_side2base(cur_kmer_info.ext_map.right);
        else
            exts.right = inv_base(ext_map_side2base(cur_kmer_info.ext_map.left));

        contig->append_base(get_base(cur_kmer, k - 1));

        if (idx == SUBCONTIG_LEN) {
            memcpy(subcontig, cur_kmer, kmer_size(k));
            idx = k;
        }
    }
    contig->exts.right = exts.right;
    return next_id;
    //while (1) {
    //    if (get_next_kmer(revcmp_found, cur_kmer_info))
    //    contig_map_type_t::iterator it = contig_map->map.find(cur_kmer);
    //    if (it == contig_map->map.end()) {
    //        revcmp_kmer(revcmp, cur_kmer, k);
    //        it = contig_map->map.find(revcmp);
    //        if (it == contig_map->map.end())
    //            break;


    //    set_base(subcontig, idx, ext_map_side2base(kmer_info.ext_map.right));
    //    idx++;
    //    for_base_in_kmer_from(b, subcontig, k, idx - k) {
    //        set_base(tmp_kmer, b_i_, b);
    //    } end_for;
    //    contig->append_kmer(tmp_kmer, k);


    //    kmer_info = it->second;

    //    if (!can_use_in_contig(kmer_info))
    //        return;
    //    /* If not reflexive */
    //    if (ext_map_side2base(kmer_info.ext_map.left) !=
    //            get_base(subcontig, idx - (k + 1)))
    //        return;

    //    if (kmer_info.contig_id >= 0) {
    //        contig->next_id = kmer_info.contig_id;
    //        return;
    //    } else {
    //        kmer_info.contig_id = contig->id;
    //    }

    //    if (idx == SUBCONTIG_LEN) {
    //        memcpy(subcontig, tmp_kmer, kmer_size(k));
    //        idx = k;
    //    }
    //}
}
#endif

void KmerCountStore::build_contigs(ContigStore& contig_store)
{
    for (contig_map_type_t::iterator it = contig_map->map.begin();
            it != contig_map->map.end();
            it++) {

        kmer_t kmer = it->first;
        kmer_info_t& kmer_info = it->second;

        if (!can_use_in_contig(kmer_info))
            continue;
        if (kmer_info.contig_found)
            continue;

        Contig* contig = new Contig(kmer);
        contig->left_ext = ext_map_side_to_base(kmer_info.ext_map.left);
        walk(contig, ext_map_side_to_base(kmer_info.ext_map.right));
        contig->revcmp();
        walk(contig, contig->right_ext);
        contig_store.add_contig(contig);
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
