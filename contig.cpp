#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "contig.h"
#include "utils.h"

using namespace std;

Contig::Contig(int32_t id)
    : subcontigs(), len(0), id(id), next_id(-1)
{
    subcontigs.push_back((kmer_t) malloc(kmer_size(SUBCONTIG_LEN)));
}

void Contig::append_kmer(kmer_t kmer, k_t len)
{
    base b;

    if (this->len == 0) {
        for_base_in_kmer(b, kmer, len) {
            set_base(subcontigs[0], b_i_, b);
        } end_for;
        this->len = len;
        return;
    }

    size_t beg_idx = (len - 1) > this->len ? 0 : this->len - (len - 1);

    uint32_t sc_idx = 0;        /* subcontig index */
    for ( ; beg_idx > SUBCONTIG_LEN; beg_idx -= SUBCONTIG_LEN)
        sc_idx++;

    /* FIXME - Get rid of these sanity checks for production code */
    if (len <= SUBCONTIG_LEN - beg_idx) {
        if (cmp_kmer(subcontigs[sc_idx], kmer, len - 1, beg_idx)) {
            char kmer_str[kmer_size(len)];
            kmer2str(kmer_str, kmer, len);
            panic("kmer (%s) didn't match when being appended to contig\n", kmer_str);
        }

        set_base(subcontigs[sc_idx], beg_idx + len - 1, get_base(kmer, len - 1));
    } else {
        size_t sublen = SUBCONTIG_LEN - beg_idx;
        size_t leftover = len - sublen;

        if (leftover == 1)
            subcontigs.push_back((kmer_t) malloc(kmer_size(SUBCONTIG_LEN)));

        if (cmp_kmer(subcontigs[sc_idx], kmer, sublen, beg_idx)) {
            char kmer_str[kmer_size(len)];
            kmer2str(kmer_str, kmer, len);
            panic("kmer (%s) didn't match when being appended to contig\n", kmer_str);
        }

        if (cmp_kmer(subcontigs[sc_idx+1], kmer, leftover - 1, 0, sublen)) {
            char kmer_str[kmer_size(len)];
            kmer2str(kmer_str, kmer, len);
            panic("kmer (%s) didn't match when being appended to contig\n", kmer_str);
        }

        set_base(subcontigs[sc_idx+1], leftover - 1, get_base(kmer, len - 1));
    }

    this->len++;
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
