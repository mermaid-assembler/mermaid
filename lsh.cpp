#include "lsh.h"
#include "config.h"
#include "base.h"

typedef uint32_t shingle_t;
static size_t SHINGLE_LEN = 16; /* This can be at most 16 */

void get_shingles(shingle_t* shingles, kmer_t kmer, size_t num_shingles)
{
    /* Initialize first shingle */
    shingles[0] = 0;
    base b;
    for_base_in_kmer(b, kmer, SHINGLE_LEN) {
        shingles[0] <<= BITS_PER_BASE;
        shingles[0] |= b;
    } end_for;

    for (size_t i = 1; i < num_shingles; i++) {
        uint32_t mask = ((1 << SHINGLE_LEN * BITS_PER_BASE) - 1) & ~3;
        uint32_t shingle = (shingles[i-1] << BITS_PER_BASE) & mask;
        shingles[i] = shingle | get_base(kmer, SHINGLE_LEN + i - 1);
    }
}

/* http://nlp.stanford.edu/IR-book/html/htmledition/near-duplicates-and-shingling-1.html */
/* Shingle + Hamming distance */
size_t lsh(kmer_t kmer, k_t k)
{
    size_t num_shingles = k - SHINGLE_LEN + 1;
    shingle_t shingles[num_shingles];
    get_shingles(shingles, kmer, num_shingles);

    uint32_t min = 0-1;
    for (size_t j = 0; j < num_shingles; j++) {
        shingle_t shingle = shingles[j];
        uint32_t idx = knuth_hash(shingle);
        if (idx < min) {
            min = idx;
        }
    }

    return min;
}

/* http://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key */
uint32_t knuth_hash(uint32_t x)
{
    return x * 2654435761;
}
