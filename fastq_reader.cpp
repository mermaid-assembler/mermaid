#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <cstring>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "fastq_reader.h"
#include "kmer.h"

using namespace std;
namespace fs = boost::filesystem;

FastQReader::FastQReader(vector<string>& fnames, k_t k)
    : k(k), paths(), curr_path(), curr_file(), so_far(), max_byte(~1), read_col(0)
{
    for (vector<string>::iterator it = fnames.begin(); it < fnames.end(); it++) {
        paths.push_back(fs::path(*it));
    }
}

uintmax_t FastQReader::total_bytes(void)
{
    uintmax_t total = 0;

    for (vector<fs::path>::iterator it = paths.begin(); it < paths.end(); it++) {
        total += fs::file_size(*it);
    }

    return total;
}

void FastQReader::set_max_byte(uintmax_t max_byte)
{
    this->max_byte = max_byte;
}

void FastQReader::seek(uintmax_t offset)
{
    for (curr_path = paths.begin(); curr_path < paths.end(); curr_path++) {
        if ((so_far + fs::file_size(*curr_path)) >= offset) {
            curr_file.close();
            curr_file.open(*curr_path);
            break;
        }
        so_far += fs::file_size(*curr_path);
    }

    /**
     * FIXME - We assume '@' will only appear at the beginning of a read. Check
     * to see that quality scores don't clash.
     */
    curr_file.seekg(offset - so_far, ios::beg);
    seek_to_next_read();
}

void FastQReader::seek_to_next_read()
{
    curr_file.ignore(MAX_LINE_LEN, '@');
    if (curr_file.eof()) {
        next_file();
    } else {
        curr_file.putback('@');
    }
}

uintmax_t FastQReader::tell(void)
{
    return so_far + curr_file.tellg();
}

bool FastQReader::next_file()
{
    curr_file.seekg(0, ios::end);
    so_far += curr_file.tellg();
    curr_path++;
    curr_file.close();
    if (curr_path == paths.end())
    {
        return false;
    }
    curr_file.open(*curr_path);
    return true;
}

bool FastQReader::read_next(qekmer_t* qekmer)
{
    while (true) {
        bool skip_kmer = false;

        if (curr_path == paths.end()) {
            return false;
        } else if (!read_col && tell() >= max_byte) {
            return false;
        }

        if (read_col ==  0) {
            //char read_buf[MAX_LINE_LEN];
            char quals_buf[MAX_LINE_LEN];

            /* Ignore the header line. */
            curr_file.ignore(MAX_LINE_LEN, '\n');
            /* Get the actual read. */
            curr_file.getline(read, MAX_LINE_LEN);
            /* Ignore +. */
            curr_file.ignore(MAX_LINE_LEN, '\n');
            /* Getthe  quality scores. */
            curr_file.getline(quals_buf, MAX_LINE_LEN);
            read_len = strlen(read);

            for (size_t i = 0; i < read_len; i++) {
                quals[i] = quals_buf[i] - ILLUMINA_QUAL_OFFSET;
            }
        }

        for (size_t i = 0; i < k; i++) {
            char char_base = read[read_col + i];
            if ( char_base == 'N') {
                read_col = read_col + i + 1;
                if (read_col > read_len - k) {
                    seek_to_next_read();
                    read_col = 0;
                }
                skip_kmer = true;
                break;
            }
            base b = char2base(char_base);
            set_base(qekmer->kmer, i, b);
        }
        if (skip_kmer) {
            continue;
        }

        if (read_col == 0) {
            qekmer->exts.left = BASE::N;
            qekmer->lqual = 0;
        } else {
            qekmer->exts.left = char2base(read[read_col - 1]);
            qekmer->lqual = quals[read_col - 1];
        }

        if (read_col + k == read_len) {
            qekmer->exts.right = BASE::N;
            qekmer->rqual = 0;
        } else {
            qekmer->exts.right = char2base(read[read_col + k]);
            qekmer->rqual = quals[read_col + k];
        }

        read_col++;

        if (read_col > read_len - k) {
            seek_to_next_read();
            read_col = 0;
        }

        return true;
    }
}
