#ifndef _LSH_H_
#define _LSH_H_

#include "kmer.h"

size_t lsh(kmer_t kmer, k_t k);
uint32_t knuth_hash(uint32_t x);

#endif /* _LSH_H_ */
