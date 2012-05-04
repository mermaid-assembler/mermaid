#ifndef _BASE_H_
#define _BASE_H_

#include <boost/cstdint.hpp>

#include "exceptions.h"

namespace BASE {
    //enum base { A = 0, C, G, T, N };
    enum base { A = 0, C, G, T };
    const uint8_t NUM_BASES = 4;

    base char2base(char c);
    char base2char(base b);
    base inv_base(base b);
    char inv_base(char c);
    bool valid_base(base b);
    bool valid_base(char c);
}

#endif /* _BASE_H_ */
