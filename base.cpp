#include "base.h"
#include "utils.h"

namespace BASE {
    base char2base(char c)
    {
        switch (c) {
            case 'A': return A;
            case 'C': return C;
            case 'G': return G;
            case 'T': return T;
            //case 'N': return N;
            case 'N': return A;
            default: panic("Invalid %s: %u\n", __func__, c);
        }
    }

    char base2char(base b)
    {
        switch(b) {
            case A: return 'A';
            case C: return 'C';
            case G: return 'G';
            case T: return 'T';
            //case N: return 'N';
            default: panic("Invalid %s: %u\n", __func__, b);
        }
    }

    base inv_base(base b)
    {
        switch (b) {
            case A: return T;
            case C: return G;
            case G: return C;
            case T: return A;
            //case N: return N;
            default: panic("Invalid %s: %u\n", __func__, b);
        }
    }

    char inv_base(char c)
    {
        switch (c) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            //case 'N': return 'N';
            default: panic("Invalid %s: %u\n", __func__, c);
        }
    }

    bool validate_base(base b)
    {
        switch (b) {
            case A:
            case C:
            case G:
            case T:
            //case N:
            //    return true;
            default:
                return false;
        }
    }
}
