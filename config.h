#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <boost/cstdint.hpp>

typedef uint32_t k_t;
typedef uint8_t qual_t;
typedef uint16_t count_t;   /* FIXME: We don't know if this is large enough for
                               the kmer counts we may see in the human dataset.
                             */


const qual_t Q_MIN = 19;   /* Minimum quality threshold. */
const count_t D_MIN = 10;   /* Minimum number of high-quality extensions needed. */

const k_t K = 41;

const unsigned int PREFIX_ALLOC_LENGTH = 5;

const size_t MIN_CONTIG_LEN = 100;
const size_t FASTA_TEXTWIDTH = 50;

#endif /* _CONFIG_H_ */
