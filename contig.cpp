#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "contig.h"
#include "utils.h"

using namespace std;

int32_t Contig::id_generator = 0;
k_t Contig::k = K;      /* TODO - Set this dynamically. */

Contig::Contig()
    : exts(), len(0), id(id_generator++), next_id(-1), subcontigs()
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

bool Contig::join_contig(Contig* next_contig)
{
    size_t subidx = len / SUBCONTIG_LEN;
    size_t sublen = len % SUBCONTIG_LEN;
    if (sublen == 0)
        subcontigs.push_back((kmer_t) malloc(kmer_size(SUBCONTIG_LEN)));

    ssize_t left_ext_len = sublen - k;
    ssize_t left_ext_idx = left_ext_len < 0 ? subidx - 1 : subidx;
    if (left_ext_len < 0) left_ext_len += SUBCONTIG_LEN;
    if (get_base(subcontigs[left_ext_idx], left_ext_len) !=
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

    return true;
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
