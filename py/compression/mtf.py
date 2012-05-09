# http://wj32.wordpress.com/2008/04/10/move-to-front-transform-in-python/
NUM_BASES = 6

def base_ord(b):
    if b == 'A': return 0
    elif b == 'C': return 1
    elif b == 'G': return 2
    elif b == 'T': return 3
    elif b == 'N': return 4
    elif b == '$': return 5
    else: raise Exception('base_ord fucked up')

def base_chr(n):
    if n == 0: return 'A'
    elif n == 1: return 'C'
    elif n == 2: return 'G'
    elif n == 3: return 'T'
    elif n == 4: return 'N'
    elif n == 5: return '$'
    else: raise Exception('base_chr fucked up')

def mtf(string):
    output=[]
    order = range(NUM_BASES)
    for c in string:
        i = order.index(base_ord(c))
        output.append(i)
        order.pop(i)
        order.insert(0, base_ord(c))
    return output

def mtf_n(string):
    '''Don't push to front if the base is N, not worth it. Difference seems to be negligible'''
    output=[]
    order = range(NUM_BASES)
    for c in string:
        i = order.index(base_ord(c))
        output.append(i)
        if i != 'N':
            order.pop(i)
            order.insert(0, base_ord(c))
    return output

def imtf(input):
    output=[]
    order = range(NUM_BASES)
    for c in input:
        output.append(base_chr(order[c]))
        order.pop(c)
        order.insert(0, base_ord(output[-1]))
    return output
