#!/usr/bin/env python

import sys

def invert_base(c):
    if c == 'A':
        return 'T'
    elif c == 'C':
        return 'G'
    elif c == 'G':
        return 'C'
    elif c == 'T':
        return 'A'
    elif c == 'N':
        return 'N'
    elif c == 'F':
        return 'F'
    elif c == 'X':
        return 'X'
    else:
        raise Exception("No reverse found for %s" % c)

def revcmp(s):
    return ''.join([invert_base(c) for c in s[::-1]])

def main():
    print(revcmp(sys.argv[1]))

if __name__ == "__main__":
    main()
