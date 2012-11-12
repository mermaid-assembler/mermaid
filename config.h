#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <boost/cstdint.hpp>
#include <string>

typedef uint32_t k_t;
typedef uint8_t qual_t;
typedef uint16_t count_t;   /* FIXME: We don't know if this is large enough for
                               the kmer counts we may see in the human dataset.
                             */

namespace Config {
    void load_config (const std::string &config_path);
    extern qual_t Q_MIN;   /* Minimum quality threshold. */
    extern count_t D_MIN;   /* Minimum number of high-quality extensions needed. */

    extern k_t K;

    extern size_t MIN_CONTIG_LEN;
    extern size_t FASTA_TEXTWIDTH;
}
#endif /* _CONFIG_H_ */
