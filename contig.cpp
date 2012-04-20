#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "contig.h"
#include "utils.h"

using namespace std;

uint32_t Contig::id_generator = 0;
k_t Contig::k = K;      /* TODO - Set this dynamically. */
size_t Contig::seed = 0;

Contig::Contig()
    : exts(), len(0), id(id_generator++), hash(), extended_hash(), revcmp_hash(),
    subcontigs()
{
    subcontigs.push_back((kmer_t) malloc(kmer_size(SUBCONTIG_LEN)));
}

Contig::~Contig()
{
    while (!subcontigs.empty()) {
        delete subcontigs.back();
        subcontigs.pop_back();
    }
}

void Contig::append_base(base b)
{
    size_t subidx = len / SUBCONTIG_LEN;
    size_t sublen = len % SUBCONTIG_LEN;
    if (sublen == 0)
        subcontigs.push_back((kmer_t) malloc(kmer_size(SUBCONTIG_LEN)));
    set_base(subcontigs[subidx], sublen, b);
    len++;
}

void Contig::append_first_kmer(kmer_t kmer)
{
    base b;
    for_base_in_kmer(b, kmer, k) {
        set_base(subcontigs[0], b_i_, b);
    } end_for;
    len = k;
}

#if 0
void Contig::append_kmer(kmer_t kmer)
{
    uint32_t subcontig_idx = get_next_free_subcontig_idx();
    size_t sublen = len - subcontig_idx * SUBCONTIG_LEN;
    size_t still_left = k;

    while (still_left > 0) {
        size_t left_in_subcontig = SUBCONTIG_LEN - sublen;
        size_t copy_len = min(still_left, left_in_subcontig);
        for_base_in_kmer(b, kmer, k
    base b;

    if (len == 0) {
        for_base_in_kmer(b, kmer, k) {
            set_base(subcontigs[0], b_i_, b);
        } end_for;
        len = k;
        return;
    }

    size_t beg_idx = (k - 1) > len ? 0 : len - (k - 1);

    uint32_t sc_idx = 0;        /* subcontig index */
    for ( ; beg_idx > SUBCONTIG_LEN; beg_idx -= SUBCONTIG_LEN)
        sc_idx++;

    /* FIXME - Get rid of these sanity checks for production code */
    if (k <= SUBCONTIG_LEN - beg_idx) {
        if (cmp_kmer(subcontigs[sc_idx], kmer, k - 1, beg_idx)) {
            char kmer_str[kmer_size(k)];
            kmer2str(kmer_str, kmer, k);
            panic("kmer (%s) didn't match when being appended to contig\n", kmer_str);
        }

        set_base(subcontigs[sc_idx], beg_idx + k - 1, get_base(kmer, k - 1));
    } else {
        size_t sublen = SUBCONTIG_LEN - beg_idx;
        size_t leftover = k - sublen;

        if (leftover == 1)
            subcontigs.push_back((kmer_t) malloc(kmer_size(SUBCONTIG_LEN)));

        if (cmp_kmer(subcontigs[sc_idx], kmer, sublen, beg_idx)) {
            char kmer_str[kmer_size(k)];
            kmer2str(kmer_str, kmer, k);
            panic("kmer (%s) didn't match when being appended to contig\n", kmer_str);
        }

        if (cmp_kmer(subcontigs[sc_idx+1], kmer, leftover - 1, 0, sublen)) {
            char kmer_str[kmer_size(k)];
            kmer2str(kmer_str, kmer, k);
            panic("kmer (%s) didn't match when being appended to contig\n", kmer_str);
        }

        set_base(subcontigs[sc_idx+1], leftover - 1, get_base(kmer, k - 1));
    }

    len++;
}
#endif

bool Contig::can_join_contig(Contig* next_contig)
{
    kmer_a extended_kmer[kmer_size(k)];
    size_t subidx = len / SUBCONTIG_LEN;
    size_t sublen = len % SUBCONTIG_LEN;

    ssize_t left_ext_sublen = sublen - k;
    ssize_t left_ext_subidx = left_ext_sublen < 0 ? subidx - 1 : subidx;
    if (left_ext_sublen < 0) left_ext_sublen += SUBCONTIG_LEN;

    if (get_base(subcontigs[left_ext_subidx], left_ext_sublen) !=
            (base) next_contig->exts.left)
        return false;

    get_extended_kmer(extended_kmer);
    if (!cmp_kmer(extended_kmer, next_contig->subcontigs[0], k))
        return true;
    else
        return false;
}

bool Contig::join_contig(Contig* next_contig)
{
    size_t subidx = len / SUBCONTIG_LEN;
    size_t sublen = len % SUBCONTIG_LEN;

    /* TODO - We can probably get rid of this chunk of code. */
    ssize_t left_ext_sublen = sublen - k;
    ssize_t left_ext_subidx = left_ext_sublen < 0 ? subidx - 1 : subidx;
    if (left_ext_sublen < 0) left_ext_sublen += SUBCONTIG_LEN;
    if (get_base(subcontigs[left_ext_subidx], left_ext_sublen) !=
            (base) next_contig->exts.left)
        return false;

    size_t next_idx = k - 1;
    size_t next_sublen = next_idx;
    size_t next_subidx = 0;
    while (next_idx < next_contig->len) {
        if (sublen == 0)
            subcontigs.push_back((kmer_t) malloc(kmer_size(SUBCONTIG_LEN)));

        base b = get_base(next_contig->subcontigs[next_subidx], next_sublen);
        set_base(subcontigs[subidx], sublen, b);

        len++;
        sublen++;
        if (sublen == SUBCONTIG_LEN) {
            subidx++;
            sublen = 0;
        }
        next_idx++;
        next_sublen++;
        if (next_sublen == SUBCONTIG_LEN) {
            next_subidx++;
            next_sublen = 0;
        }
    }

    extended_hash = next_contig->extended_hash;
    revcmp_hash = next_contig->revcmp_hash;

    return true;
}

void Contig::revcmp(void)
{
    kmer_a extended_kmer[kmer_size(k)];
    kmer_a revcmp[kmer_size(len)];
    base b;
    size_t subidx;
    size_t sublen;

    for (subidx = 0; subidx < len / SUBCONTIG_LEN; subidx++) {
        for_base_in_kmer(b, subcontigs[subidx], SUBCONTIG_LEN) {
            set_base(revcmp, (len - 1) - (subidx * SUBCONTIG_LEN + b_i_), inv_base(b));
        } end_for;
    }

    sublen = len % SUBCONTIG_LEN;
    for_base_in_kmer(b, subcontigs[subidx], sublen) {
        set_base(revcmp, (len - 1) - (subidx * SUBCONTIG_LEN + b_i_), inv_base(b));
    } end_for;

    for_base_in_kmer(b, revcmp, len) {
        subidx = b_i_ / SUBCONTIG_LEN;
        sublen = b_i_ % SUBCONTIG_LEN;
        set_base(subcontigs[subidx], sublen, b);
    } end_for;

    get_extended_kmer(extended_kmer);
    extended_hash = kmer_hash(seed, extended_kmer, k);
    swap(hash, revcmp_hash);
}


void Contig::generate_hashes(void)
{
    hash = kmer_hash(seed, subcontigs[0], k);

    kmer_a extended_kmer[kmer_size(k)];
    kmer_a revcmp[kmer_size(k)];

    for (size_t i = 0; i < k; i++) {
        size_t subidx = (len - k + i) / SUBCONTIG_LEN;
        size_t sublen = (len - k + i) % SUBCONTIG_LEN;

        base b = get_base(subcontigs[subidx], sublen);
        set_base(revcmp, (k - 1) - i, inv_base(b));
        if (i > 0)
            set_base(extended_kmer, i - 1, b);
    }
    set_base(extended_kmer, k - 1, (base) exts.right);

    extended_hash = kmer_hash(seed, extended_kmer, k);
    revcmp_hash = kmer_hash(seed, revcmp, k);
}

void Contig::fprint(FILE* outfile)
{
    size_t still_left = len;
    uint32_t sc_idx = 0;
    while (still_left > 0) {
        size_t to_print = min(still_left, SUBCONTIG_LEN);
        fprint_kmer(outfile, subcontigs[sc_idx], to_print);
        still_left -= to_print;
        sc_idx++;
    }
}

void Contig::fprintln(FILE* outfile)
{
    fprint(outfile);
    fprintf(outfile, "\n");
}

void Contig::get_extended_kmer(kmer_t extended_kmer)
{
    for (size_t i = 0; i < k - 1; i++) {
        size_t subidx = (len - (k - 1) + i) / SUBCONTIG_LEN;
        size_t sublen = (len - (k - 1) + i) % SUBCONTIG_LEN;

        base b = get_base(subcontigs[subidx], sublen);
        set_base(extended_kmer, i, b);
    }
    set_base(extended_kmer, k - 1, (base) exts.right);
}
