#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "kmer.h"
#include "contig.h"
#include "utils.h"

using namespace std;

uint32_t Contig::id_generator = 0;
k_t Contig::k = K;      /* TODO - Set this dynamically. */
size_t Contig::seed = 0;

Contig::Contig(kmer_t init_kmer)
    : id(id_generator++), s(), left_ext(), right_ext()
{
    char kmer_str[kmer_size(k)];
    kmer2str(kmer_str, init_kmer, k);
    s = kmer_str;
}

bool Contig::check_next_left_ext(base next_left_ext)
{
    return char2base(s[s.size() - k]) == next_left_ext;
}

void Contig::revcmp(void)
{
    reverse(s.begin(), s.end());
    for (string::iterator it = s.begin(); it < s.end(); it++)
        *it = inv_base(*it);

    swap(left_ext, right_ext);
    left_ext = inv_base(left_ext);
    right_ext = inv_base(right_ext);
}

void Contig::fprint(FILE* outfile)
{
    fprintf(outfile, "%s", s.c_str());
}

void Contig::fprintln(FILE* outfile)
{
    fprint(outfile);
    fprintf(outfile, "\n");
}
