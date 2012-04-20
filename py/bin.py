#!/usr/bin/env python

import sys
import random

import kmerutils

NUM_KMERS = 10000
K = 41
NUM_BINS = 16
random.seed(1073847)
SEEDS = [random.randint(0, sys.maxint) for i in range(4)]

def my_hash(kmer):
    def naive_hash(kmer):
        return {
                'AA':  0,
                'AC':  1,
                'AG':  2,
                'AT':  3,
                'CA':  4,
                'CC':  5,
                'CG':  6,
                'CT':  7,
                'GA':  8,
                'GC':  9,
                'GG': 10,
                'GT': 11,
                'TA': 12,
                'TC': 13,
                'TG': 14,
                'TT': 15,
                }[kmer[:2].replace('N', 'A')]
    def naive_hash2(kmer):
        return hash(kmer) % NUM_BINS
    def lsh(kmer):
        # Kind of bullshitting my way through this from info here:
        # http://nlp.stanford.edu/IR-book/html/htmledition/near-duplicates-and-shingling-1.html
        def permute(x, seed):
            random.seed(seed)
            return random.randint(0, x)
        # Params to do some search exploration over:
        # * shingle length
        # * number of seeds to use
        # * could we select the shingles in a strided pattern?
        shingle_len = 6 # because 2^6 == 64, 1 bit per shingle
        shingles = set()
        for i in range(len(kmer)-shingle_len+1):
            shingles.add(kmerutils.int_value(kmer[i:i+6]))
        results = []
        for seed in SEEDS:
            min_permuted = min([permute(x, seed) for x in shingles])
            results.append(min_permuted)
        value = 0
        for result in results:
            value *= 2
            value += result % 2
        return value
    return lsh(kmer)

def main():
    if len(sys.argv) != 2:
        print('Usage: %s fastq-file' % sys.argv[0])
        sys.exit(1)

    k = K
    reader = kmerutils.FastQReader(sys.argv[1])
    bins = []
    for i in range(NUM_BINS):
        bins.append([])

    counts = [0]*NUM_BINS
    for i in range(NUM_KMERS):
        kmer = reader.next_kmer(k+1)
        # Skip Ns
        while kmer.find('N') != -1:
            kmer = reader.next_kmer(k+1)
        bin_id = my_hash(kmer[:-1])
        # Just count forward extensions WLOG
        if bin_id == my_hash(kmer[1:]):
            counts[bin_id] += 1
        bins[bin_id].append(kmer[1:-1])

    print("bin sizes:           %s" % " ".join(["%6d" % len(bin) for bin in bins]))

    print("intrabin extensions: %s" % " ".join(["%6d" % c for c in counts]))

    locality = [(100.0*x/y if y > 0 else 0) for x,y in zip(counts, [len(bin) for bin in bins])]
    print("locality:             %s" % " ".join(["%5.2f%%" % x for x in locality]))

    print("average locality: %5.2f" % (100.0 * sum(counts)/sum([len(bin) for bin in bins])))


if __name__ == '__main__':
    main()
