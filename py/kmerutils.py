#!/usr/bin/env python

import sys

K = 41
ILLUMINA_QUAL_OFFSET = 64
Q_MIN = 20

NUM_BASES = 6

def base_ord(b):
    if b == 'A': return 0
    elif b == 'C': return 1
    elif b == 'G': return 2
    elif b == 'T': return 3
    elif b == 'N': return 4
    elif b == '$': return 5
    else: raise Exception("base_ord got recognized base: '%s'" % b)

def base_chr(n):
    if n == 0: return 'A'
    elif n == 1: return 'C'
    elif n == 2: return 'G'
    elif n == 3: return 'T'
    elif n == 4: return 'N'
    elif n == 5: return '$'
    else: raise Exception('base_chr fucked up')

def inv_base(c):
    if c == 'A': return 'T'
    elif c == 'C': return 'G'
    elif c == 'G': return 'C'
    elif c == 'T': return 'A'
    elif c == 'N': return 'N'
    elif c == 'B': return 'B'
    else:
        raise Exception("No reverse found for %s" % c)

def int_value(kmer):
    '''Returns a unique integer representing the kmer.'''
    num = 0
    for base in kmer:
        num *= 4
        num += base_ord(base)
    return num

def revcmp_kmer(kmer):
    return ''.join([inv_base(b) for b in kmer[::-1]])

def canonize(kmer):
    revcmp = revcmp_kmer(kmer)
    if revcmp < kmer:
        return revcmp
    else:
        return kmer

def score_kmer(kmer):
    score = 0
    for base in kmer:
        score = score * 4 + base_ord(base)
    return score

def rotate_kmer(kmer, num, k = K):
    rot = kmer
    for i in range(num):
        rot = rot[1:] + rot[k - 1]
    return rot

def canonize_cycle(kmer, k = K):
    if kmer[:k - 1] == kmer[-(k - 1):]:
        # It's a cycle
        score_map = {}
        for i in range(len(kmer)):
            rot = rotate_kmer(kmer, i)
            score_map[score_kmer(rot)] = rot
        revcmp = revcmp_kmer(kmer)
        for i in range(len(revcmp)):
            rot = rotate_kmer(revcmp, i)
            score_map[score_kmer(rot)] = rot
        return score_map[sorted(score_map.keys())[0]]
    else:
        return canonize(kmer)

class FastAReader(object):
    def __init__(self, fname):
        self.f = open(fname)
        self.f.readline()

    def next_contig(self):
        contig = ''
        while True:
            read = self.f.readline().strip()
            if read == '' or read[0] == '>':
                return contig
            contig += read

    def contigs(self):
        contig = ''
        for line in self.f.readlines():
            read = line.strip()
            if read == '':
                continue
            elif read[0] == '>':
                yield contig
                contig = ''
            else:
                contig += read
        if contig != '':
            yield contig

class FastQReader(object):
    def __init__(self, fname):
        self.f = open(fname, 'r')
        self.read = ''
        self.read_idx = 0

    def next_read(self):
        self.f.readline()
        read = self.f.readline().strip()
        self.f.readline()
        quals = self.f.readline().strip()
        return read, quals

    def next_kmer(self, k):
        if self.read_idx + k >= len(self.read):
            self.read, quals = self.next_read()
            self.read_idx = 0
        
        if self.read == '':
            return ''

        kmer = self.read[self.read_idx:self.read_idx + k]
        self.read_idx += 1

        return kmer

# Returns a list of (read, quals) tuples which contain the kmer or the reverse
# complement of the kmer. Actually, this is implemented as a generator
def containing_reads(kmer, reader):
    revcmp = revcmp_kmer(kmer)
    while True:
        read, quals = reader.next_read()
        if read == '':
            break
        if kmer in read or revcmp in read:
            yield (read, quals)

def main():
    if len(sys.argv) != 3:
        print("Usage: %s kmer filename" % sys.argv[0])
        sys.exit(2)
    kmer = sys.argv[1]
    fname = sys.argv[2]
    reader = FastQReader(fname)

    for read, quals in containing_reads(kmer, reader):
        print '@HWI-EAS306:7:1:1003:1817#NNNNNN/1'
        print read
        print '+'
        print quals

if __name__ == '__main__':
    main()
