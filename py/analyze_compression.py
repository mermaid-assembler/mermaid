import random
import sys

import compression.bwt as bwt
import compression.mtf as mtf
import compression.encoding as enc
import kmerutils

NUM_KMERS = 10000
#NUM_KMERS = 1
K = 41

def get_random_kmers(reader, num, k):
    kmers = []
    total_kmers = num*100
    for i in range(total_kmers):
        next_kmer = reader.next_kmer(k)
        # skip kmers with Ns
        #while next_kmer.find('N') != -1:
        #    next_kmer = reader.next_kmer(k)
        kmers.append(next_kmer)
    for i in range(total_kmers):
        idx = random.randint(0,i)
        idx2 = total_kmers - i - 1
        tmp = kmers[idx2]
        kmers[idx2] = kmers[idx]
        kmers[idx] = tmp
    return kmers[:num]

def convert_to_rle(compressed):
    rle = [(compressed[0], 1)]
    i = 1
    while i < len(compressed):
        if compressed[i] == rle[-1][0]:
            rle[-1] = (rle[-1][0], rle[-1][1] + 1)
        else:
            rle.append((compressed[i], 1))
        i += 1
    return rle

    #def printable_rle(rle):
    #    s = ''
    #    for elem, num in rle:
    #        s += '%d:%d ' % (elem, num)
    #    return s.rstrip()

    #rle = convert_to_rle(mtf.mtf(bwt.bwt(kmer)))
    #rle = convert_to_rle([kmerutils.base_ord(b) for b in bwt.bwt(kmer)])
    #print '%s\t%s' % (kmer, printable_rle(rle))

def main():
    if len(sys.argv) != 2:
        print('Usage: %s datafile' % sys.argv[0])
        sys.exit(2)

    reader = kmerutils.FastQReader(sys.argv[1])
    kmers = get_random_kmers(reader, NUM_KMERS, K)
    runs = {}
    for kmer in kmers:
        kmer = kmer.strip()
        compressed = [kmerutils.base_ord(b) for b in bwt.bwt(kmer)]
        #compressed = mtf.mtf(bwt.bwt(kmer))
        for elem, num in convert_to_rle(compressed):
            runs[num] = runs.get(num, 0) + 1

    total_bases = NUM_KMERS * K
    lengths = sorted(runs.keys())
    for i in lengths:
        print '%d: %d (%f)' % (i, runs[i], 100 * i * float(runs[i]) / total_bases)
    savings = sum([(i - 2) * runs[i] for i in lengths if i > 2])
    print 'total bases: %d savings: %d percent: %f' % \
            (total_bases, savings, 100 * float(savings) / total_bases)

    # Size in bytes assuming 2 bits per base and blocks of 8 bits
    bits_per_kmer = 2 * K
    bytes_per_kmer = (bits_per_kmer+1)/8 - 1
    total_size = bytes_per_kmer * NUM_KMERS
    total_bits = bits_per_kmer * NUM_KMERS
    print("Naive encoding size (bytes, bits): %d, %d" % (total_size, total_bits))

    bits_per_kmer = 3 * K
    bytes_per_kmer = (bits_per_kmer+1)/8 - 1
    total_size = bytes_per_kmer * NUM_KMERS
    total_bits = bits_per_kmer * NUM_KMERS
    print("Naive encoding size (3 bits) (bytes, bits): %d, %d" % (total_size, total_bits))
    naive_size = total_size

    size = 0
    bit_size = 0
    for kmer in kmers:
        rle_bits = enc.scheme6(kmer)
        bit_size += rle_bits
        rle_bytes = (rle_bits + 1)/8 - 1
        size += rle_bytes
    print("encoding size (bytes, bits): %d, %d" % (size, bit_size))
    print("%f%% savings" % (100 * (naive_size-size)*1.0/naive_size))

if __name__ == '__main__':
    main()
