#ifndef _FASTQ_READER_H_
#define _FASTQ_READER_H_

#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "kmer.h"

/* FIXME - This should be more dynamic. */
const int MAX_LINE_LEN = 1024;


const int ILLUMINA_QUAL_OFFSET = 64;

class FastQReader {
public:
    FastQReader(std::vector<std::string>& fnames, k_t k);

    uintmax_t total_bytes(void);
    void seek(uintmax_t offset);
    void set_max_byte(uintmax_t max_byte);
    uintmax_t tell(void);
    bool next_file();
    bool read_next(qekmer_t* qekmer);

protected:
    void seek_to_next_read();

    k_t k;

    std::vector<boost::filesystem::path> paths;
    std::vector<boost::filesystem::path>::iterator curr_path;
    boost::filesystem::ifstream curr_file;
    uintmax_t so_far;           // At the granularity of files
    uintmax_t max_byte;

    char read[MAX_LINE_LEN];    /* In 'base' type. */
    qual_t quals[MAX_LINE_LEN];
    size_t read_col;            /* The column index within the current read that
                                   the next kmer (acutal kmer; not k+2-mer)
                                   starts from. */
    size_t read_len;
};

#endif /* _FASTQ_READER_H_ */
