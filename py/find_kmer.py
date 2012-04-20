#!/usr/bin/env python

import sys
import re
from revcmp import revcmp, invert_base

ILLUMINA_QUAL_OFFSET = 64
Q_MIN = 20

class FastQReader(object):
    def __init__(self, file):
        self.f = open(file, 'r')

    def read_next(self):
        self.f.readline()
        read = self.f.readline().strip()
        self.f.readline()
        quals = self.f.readline().strip()
        return read, quals

def find_all(read, mer):
    start = 0
    #lst = []
    while True:
        idx = read.find(mer, start)
        if idx != -1:
            #lst.append(idx)
            yield idx
            start = idx + 1
        else:
            break
    #return lst
count = 0
good_l_count = 0
good_r_count = 0

def look_for_extensions(mer, reader):
    def find(read, quals, mer, reverse=False):
        global count, good_l_count, good_r_count
        for idx in find_all(read, mer):
            count += 1
            if idx > 0:
                l = read[idx-1] 
                l_qual = ord(quals[idx-1]) - ILLUMINA_QUAL_OFFSET
            else:
                l = 'X'
                l_qual = 0

            if idx + len(mer) < len(read):
                r = read[idx + len(mer)] 
                r_qual = ord(quals[idx + len(mer)]) - ILLUMINA_QUAL_OFFSET
            else:
                r = 'X'
                r_qual = 0

            if reverse:
                if r_qual >= Q_MIN:
                    good_l_count += 1
                if l_qual >= Q_MIN:
                    good_r_count += 1
                print('%s%s\t(%d,%d)\tR' % (invert_base(r),invert_base(l),r_qual,l_qual))
            else:
                if l_qual >= Q_MIN:
                    good_l_count += 1
                if r_qual >= Q_MIN:
                    good_r_count += 1
                print('%s%s\t(%d,%d)' % (l,r,l_qual,r_qual))

    rev = revcmp(mer)
    while True:
        read, quals = reader.read_next()
        if read == "":
            break
        find(read, quals, mer)
        find(read, quals, rev, reverse=True)
    print("Total lines: %d" % count)
    print("Good left extensions: %d" % good_l_count)
    print("Good right extensions: %d" % good_r_count)

def main():
    if len(sys.argv) != 3:
        print("Usage: %s mer filename" % sys.argv[0])
        sys.exit(1)
    mer = sys.argv[1]
    file = sys.argv[2]
    reader = FastQReader(file)
    look_for_extensions(mer, reader)

if __name__ == "__main__":
    main()
