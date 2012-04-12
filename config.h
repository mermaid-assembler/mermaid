#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <boost/cstdint.hpp>

typedef uint32_t k_t;

const uint8_t Q_MIN = 19;   /* Minimum quality threshold. */
const uint8_t D_MIN = 10;   /* Minimum number of high-quality extensions needed. */

const k_t K = 41;

const unsigned int PREFIX_ALLOC_LENGTH = 5;

#endif /* _CONFIG_H_ */
