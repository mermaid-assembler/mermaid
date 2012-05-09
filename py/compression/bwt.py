#http://barnesc.blogspot.com/2007/04/burrows-wheeler-transform-in-python.html
def bwt(s):
  s = s + '$'
  return ''.join([x[-1] for x in
         sorted([s[i:] + s[:i] for i in range(len(s))])])

def ibwt(s):
  L = [''] * len(s)
  for i in range(len(s)):
    L = sorted([s[i] + L[i] for i in range(len(s))])
  return [x for x in L if x.endswith('$')][0][:-1]

# For a reduced alphabet size, we don't store the end character in the string
# Example: AB$ADC -> (ABADC, 2)
def eol_format(s):
    idx = s.find('$')
    new_s = [x for x in s]
    new_s.pop(idx)
    return (new_s, idx)

def ieol_format(s):
    string = s[0]
    string.insert(s[1], '$')
    return ''.join(string)
