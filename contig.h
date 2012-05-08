#ifndef _CONTIG_H_
#define _CONTIG_H_

#include <cstdio>
#include <vector>
#include <string>

#include <boost/cstdint.hpp>

#include "kmer.h"

class Contig {
public:
    static void set_k(k_t k) { Contig::k = k; }

    Contig();
    Contig(kmer_t init_kmer);
    Contig(Contig* contig);

    /* Checks the given base 'next_left_ext' against the s[size - k] to see if
     * they're the same. If the kmer overlaps properly, these bases should be
     * the same.
     */
    bool check_next_left_ext(base next_left_ext);

    /* Reverse complements the current contig. */
    void revcmp(void);

    /* Sets kmer to [last k-1 bases of s] + [right_ext] */
    void get_ext_kmer(kmer_t kmer);

    void fprint(FILE* outfile);
    void fprintln(FILE* outfile);
    /* Print formatted contig to file. */
    void fprint_fasta(FILE* outfile, size_t textwidth);

    /* Append an overlapping contig. */
    void append(Contig* next_contig);

    /* Checks for badness */
    void verify();

    /* Checks whether contigs contains the given kmer. */
    bool contains(const char* kmer_str);
    bool contains_kmer_or_revcmp(const char* kmer_str);

public:
    uint32_t id;
    std::string s;      /* string that represents the contig. */
    base left_ext;
    base right_ext;

protected:
    static uint32_t id_generator;
    static k_t k;
    static size_t seed;
};

#endif /* _CONTIG_H_ */
