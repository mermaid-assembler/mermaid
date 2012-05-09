import bwt
import mtf

def run_len(kmer, idx):
    i = idx
    size = 0
    c = kmer[idx]
    while i < len(kmer):
        if kmer[i] == c:
            size += 1
        else:
            break
        i += 1
    return size

# These guys all return number of bits
def scheme1(kmer):
    '''
    BWT, EOL, 2 bit alphabet
    Encoding scheme = 
    encoding := eol_pos bases
    eol_pos := 6 bits for eol, suitable for k <=~ 63
    bases := 2 bits for base, 1 bit for RLE, 1 bit for len 2-3, (1 bit for len 2-3 | 6 bits for len <= 63)
    Examples: AC == 000 010
              AAC == 00100 010
              AAAC == 00101 010
              AAAAC == 0011000000 010
              AAAAAC == 0011000001 010
    Three-tiered encoding, run len 1 => 3 bits, run len 2-3 => 5 bits, run len 4-68 => 10 bits
    '''
    trans_kmer, eol_idx = bwt.eol_format(bwt.bwt(kmer))
    size = 0
    size += 6
    i = 0
    while i < len(trans_kmer):
        run_length = run_len(trans_kmer, i)
        if run_length == 1:
            size += 3
        elif run_length in [2,3]:
            size += 5
        else:
            size += 10
        if run_length > 68:
            run_length = 68
        i += run_length
    return size

def scheme2(kmer):
    '''
    BWT, 3 bit alphabet
    Encoding scheme = 
    encoding := bases
    bases := 3 bits for base, 1 bit for RLE, 4 bits for len 2-17
    Examples: AC     == 0000 0010
              AAC    == 00010000 0010
              AAAC   == 00010001 0010
              AAAAC  == 00010010 0010
              AAAAAC == 00010011 0010
    Two-tiered encoding, run len 1 => 4 bits, run len 2-17 => 8 bits
    '''
    trans_kmer = bwt.bwt(kmer)
    size = 0
    i = 0
    while i < len(trans_kmer):
        run_length = run_len(trans_kmer, i)
        if run_length == 1:
            size += 4
        else:
            size += 8
        if run_length > 17:
            run_length = 17
        i += run_length
    return size

def scheme3(kmer):
    '''
    2 bit alphabet, BWT, MTF, shitty huffman encoding, 0 => 1, 1 => 01, 2 => 000, 3 => 001
    Add 6 bits for eol position
    '''
    def huff_size(c):
        return { 0: 1, 
                 1: 2,
                 2: 3,
                 3: 3,
                 }[c]
    trans_kmer, eol_idx = bwt.eol_format(bwt.bwt(kmer))
    trans_kmer = mtf.mtf(trans_kmer)
    size = 0
    size += 6
    i = 0
    while i < len(trans_kmer):
        size += huff_size(trans_kmer[i])
        i += 1
    return size

def scheme4(kmer):
    '''
    3 bit alphabet, BWT, MTF, shitty huffman encoding, 0 => 1, 1 => 01, 2 => 000, 3 => 0011, 4 => 0010
    Add 6 bits for eol position
    '''
    def huff_size(c):
        return { 0: 1, 
                 1: 2,
                 2: 3,
                 3: 4,
                 4: 4,
                 }[c]
    trans_kmer, eol_idx = bwt.eol_format(bwt.bwt(kmer))
    trans_kmer = mtf.mtf(trans_kmer)
    size = 0
    size += 6
    i = 0
    while i < len(trans_kmer):
        size += huff_size(trans_kmer[i])
        i += 1
    return size

def scheme5(kmer):
    '''
    3 bit alphabet, BWT, MTF, shitty huffman encoding, 0 => 1, 1 => 01, 2 => 000, 3 => 0011, 4 => 0010
    Add 6 bits for eol position
    modified MTF to not push to front when base is N
    '''
    def huff_size(c):
        return { 0: 1, 
                 1: 2,
                 2: 3,
                 3: 4,
                 4: 4,
                 }[c]
    trans_kmer, eol_idx = bwt.eol_format(bwt.bwt(kmer))
    trans_kmer = mtf.mtf_n(trans_kmer)
    size = 0
    size += 6
    i = 0
    while i < len(trans_kmer):
        size += huff_size(trans_kmer[i])
        i += 1
    return size

def scheme6(kmer):
    '''
    3 bit alphabet, MTF, shitty huffman encoding, 0 => 1, 1 => 01, 2 => 000, 3 => 0011, 4 => 0010
    modified MTF to not push to front when base is N
    '''
    def huff_size(c):
        return { 0: 1, 
                 1: 2,
                 2: 3,
                 3: 4,
                 4: 4,
                 }[c]
    trans_kmer = mtf.mtf_n(kmer)
    size = 0
    i = 0
    while i < len(trans_kmer):
        size += huff_size(trans_kmer[i])
        i += 1
    return size

