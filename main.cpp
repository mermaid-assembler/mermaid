#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <map>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/mpi.hpp>

#include "fastq_reader.h"
#include "nethub.h"
#include "config.h"
#include "hash_map.h"
#include "scalable_bloom_filter.h"

using namespace std;
namespace mpi = boost::mpi;
namespace fs  = boost::filesystem;

#define FILTER_ON 0

#define KMER_BUFFER_SIZE 1024
#define KMER_TAG 0
#define KMER_SIZE_TAG 1
#define DONE_TAG 2

static const k_t k = K;

/* FIXME - Change this initial capacity using preprocessing step. */
static const int INITIAL_CAPACITY = 100000;

int get_kmer_bin(qkmer_t* qkmer, k_t k, int world_size)
{
    base b;
    kmer_a kmer[kmer_size(k)];
    for_base_in_kmer_from(b, qkmer->kmer, k, 1) {
        set_base(kmer, b_i_, b);
    } end_for;
    size_t hash = kmer_hash(0, kmer, k);
    return hash % world_size;
}

/* Checks to see if the qkmer meets the criteria to being sent. For example,
 * this function checks to make sure that the kmer does not contain an 'N'
 * anywhere in the kmer.
 */
bool check_qkmer_qual(qkmer_t* qkmer, k_t k)
{
    base b;
    for_base_in_kmer_from(b, qkmer->kmer, k, 1) {
        if (b == BASE::N)
            return false;
    } end_for;

    if (get_base(qkmer->kmer, 0) == BASE::N &&
            get_base(qkmer->kmer, k + 1) == BASE::N)
        return false;

    if (qkmer->lqual <= Q_MIN && qkmer->rqual <= Q_MIN)
        return false;

    return true;
}

/* Sets the qkmer to the canonical qkmer. */
void canonize_qkmer(qkmer_t* qkmer, k_t k)
{
    kmer_a revcmp[kmer_size(k + 2)];
    revcmp_kmer(revcmp, qkmer->kmer, k + 2);
    if (cmp_kmer(qkmer->kmer, revcmp, k, 1, 1) > 0) {
        memcpy(qkmer->kmer, revcmp, kmer_size(k + 2));
        qual_t tmp = qkmer->lqual;
        qkmer->lqual = qkmer->rqual;
        qkmer->rqual = tmp;
    }
}

/* Replaces endings with N's to A's with 0 qual score. */
void clean_qkmer(qkmer_t* qkmer, k_t k)
{
    if (get_base(qkmer->kmer, 0) == BASE::N) {
        set_base(qkmer->kmer, 0, BASE::A);
        qkmer->lqual = 0;
    }
    if (get_base(qkmer->kmer, k+1) == BASE::N) {
        set_base(qkmer->kmer, k+1, BASE::A);
        qkmer->rqual = 0;
    }
}

FastQReader* get_reader(int argc, char* argv[], mpi::communicator& world, k_t k)
{
    vector<string> fnames;
    for (int i = 0; i < argc; ++i) {
        fnames.push_back(string(argv[i]));
    }

    FastQReader* reader = new FastQReader(fnames, k);
    uintmax_t offset = reader->total_bytes() / world.size() * world.rank();
    reader->seek(offset);
    reader->set_max_byte(offset + reader->total_bytes() / world.size());

    return reader;
}

void build_map(FastQReader* r, HashMap<kmer_t, qual_counts_t>& kmer_map, mpi::communicator& world)
{
    NetHub nethub(world, qkmer_size(k+2));

    ScalableBloomFilter kmer_filter(INITIAL_CAPACITY, 0.001, 0,
            (bloom_filter_hash_func_t) kmer_hash_K);

    bool node_done[world.size()]; // Stores the non-blocking receive requests
    for (int i = 0; i < world.size(); i++) {
        node_done[i] = false;
    }
    bool done_reading = false;
    bool all_done = false;

    qkmer_t* send_qkmer = (qkmer_t*) malloc(qkmer_size(k + 2));
    qkmer_t* recv_qkmer = (qkmer_t*) malloc(qkmer_size(k + 2));
    base b;
    kmer_a recv_kmer[kmer_size(k)];

    while (!all_done) {
        if (r->read_next(send_qkmer)) {
            if (check_qkmer_qual(send_qkmer, k)) {
                canonize_qkmer(send_qkmer, k);
                clean_qkmer(send_qkmer, k);
                int node_id = get_kmer_bin(send_qkmer, k, world.size());
                nethub.send(node_id, send_qkmer);
            }
        } else {
            if (!done_reading) {
                nethub.done();
                done_reading = true;
            }
        }

        all_done = true;
        for (int i = 0; i < world.size(); i++)
        {
            if (node_done[i]) continue;

            int status;
            while ((status = nethub.recv(i, recv_qkmer)) == 0) {
                for_base_in_kmer_from(b, recv_qkmer->kmer, k, 1) {
                    set_base(recv_kmer, b_i_, b);
                } end_for;

#if FILTER_ON
                if (!kmer_filter.check(recv_kmer)) {
                    kmer_filter.add(recv_kmer);
                } else {
#endif
                    kmer_map.try_insert(recv_kmer);
                    if (recv_qkmer->lqual > Q_MIN)
                        kmer_map.map[recv_kmer].lquals[get_base(recv_qkmer->kmer, 0)]++;
                    if (recv_qkmer->rqual > Q_MIN)
                        kmer_map.map[recv_kmer].rquals[get_base(recv_qkmer->kmer, k+1)]++;
#if FILTER_ON
                }
#endif
            }

            if (status == 1)
                node_done[i] = true;
            else
                all_done = false;
        }
    }

    free(send_qkmer);
    free(recv_qkmer);
}

static inline bool check_hq_ext(ext_t ext)
{
    /* This kind of assumes that ext_t := uint8 */
    return ext.bitmap != 0;
    /*
    return ext.lA || ext.lC || ext.lG || ext.lT ||
           ext.rA || ext.rC || ext.rG || ext.rT;
    */
}

ext_t qual_counts_2_ext(qual_counts_t* qual_counts)
{
    ext_t ext;
    ext.bitmap = 0;
    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        if (qual_counts->lquals[i] >= D_MIN) {
            ext.bitmap |= 1 << (4 + i);
        }
        if (qual_counts->rquals[i] >= D_MIN) {
            ext.bitmap |= 1 << i;
        }
    }
    return ext;
}

void gather_kmers(HashMap<kmer_t, ext_t>& all_kmer_map,
        HashMap<kmer_t, qual_counts_t>& kmer_map, mpi::communicator& world)
{
    ekmer_t* send_ekmer = (ekmer_t*) malloc(ekmer_size(k));
    ekmer_t* recv_ekmer = (ekmer_t*) malloc(ekmer_size(k));
    base b;

    NetHub nethub(world, ekmer_size(k));

    bool node_done[world.size()]; // Stores the non-blocking receive requests
    for (int i = 0; i < world.size(); i++) {
        node_done[i] = false;
    }
    bool done_sending = false;
    bool all_done = false;

    HashMap<kmer_t, qual_counts_t>::map_type_t::iterator it = kmer_map.map.begin();
    while (!all_done) {
        if (it != kmer_map.map.end()) {
            const kmer_t& kmer = it->first;
            qual_counts_t& qual_counts = it->second;
            ext_t ext = qual_counts_2_ext(&qual_counts);

            if (check_hq_ext(ext)) {
                if (world.rank() == 0) {
                    all_kmer_map.try_insert(kmer);
                    all_kmer_map.map[kmer] = ext;
                } else {
                    for_base_in_kmer(b, kmer, k) {
                        set_base(send_ekmer->kmer, b_i_, b);
                    } end_for;
                    send_ekmer->ext = ext;
                    nethub.send(0, send_ekmer);
                }
            }

            it++;
        } else {
            if (!done_sending) {
                nethub.done();
                done_sending = true;
            }
        }
        
        if (world.rank() != 0) {
            all_done = done_sending;
        } else {
            all_done = true;
            for (int i = 0; i < world.size(); i++) {
                if (node_done[i]) continue;

                int status;
                while ((status = nethub.recv(i, recv_ekmer)) == 0) {
                    all_kmer_map.try_insert(recv_ekmer->kmer);
                    all_kmer_map.map[recv_ekmer->kmer] = recv_ekmer->ext;
                }

                if (status == 1)
                    node_done[i] = true;
                else
                    all_done = false;
            }
        }
    }

    free(send_ekmer);
    free(recv_ekmer);
}

void distrib_print_ufx(FILE* outfile, HashMap<kmer_t, qual_counts_t> kmer_map)
{
    for (HashMap<kmer_t, qual_counts_t>::map_type_t::iterator it = kmer_map.map.begin();
            it != kmer_map.map.end();
            it++) {
        const kmer_t& kmer = it->first;
        qual_counts_t& qual_counts = it->second;
        ext_t ext = qual_counts_2_ext(&qual_counts);

        if (!check_hq_ext(ext))
            continue;

        char left_ext = 0;
        char right_ext = 0;

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (ext.bitmap & 1 << (BASE::NUM_BASES + i)) {
                if (left_ext)
                    left_ext = 'F';
                else 
                    left_ext = base2char((base) i);
            }
        }
        if (!left_ext)
            left_ext = 'X';

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (ext.bitmap & 1 << i) {
                if (right_ext)
                    right_ext = 'F';
                else 
                    right_ext = base2char((base) i);
            }
        }
        if (!right_ext)
            right_ext = 'X';


        fprint_kmer(outfile, kmer, k);
        fprintf(outfile, "\t%c%c\n", left_ext, right_ext);

        kmer_a revcmp[kmer_size(k)];
        revcmp_kmer(revcmp, kmer, k);
        fprint_kmer(outfile, revcmp, k);
        fprintf(outfile, "\t");
        if (right_ext == 'F' || right_ext == 'X')
            fprintf(outfile, "%c", right_ext);
        else
            fprintf(outfile, "%c", inv_base(right_ext));
        if (left_ext == 'F' || left_ext == 'X')
            fprintf(outfile, "%c", left_ext);
        else
            fprintf(outfile, "%c", inv_base(left_ext));
        fprintf(outfile, "\n");
    }

}

void print_ufxs(FILE* outfile, HashMap<kmer_t, ext_t>& all_kmer_map)
{
    for (HashMap<kmer_t, ext_t>::map_type_t::iterator it = all_kmer_map.map.begin();
            it != all_kmer_map.map.end();
            it++) {

        char left_ext = 0;
        char right_ext = 0;

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (it->second.bitmap & 1 << (BASE::NUM_BASES + i)) {
                if (left_ext)
                    left_ext = 'F';
                else 
                    left_ext = base2char((base) i);
            }
        }
        if (!left_ext)
            left_ext = 'X';

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (it->second.bitmap & 1 << i) {
                if (right_ext)
                    right_ext = 'F';
                else 
                    right_ext = base2char((base) i);
            }
        }
        if (!right_ext)
            right_ext = 'X';


        fprint_kmer(outfile, it->first, k);
        fprintf(outfile, "\t%c%c\n", left_ext, right_ext);

        kmer_a revcmp[kmer_size(k)];
        revcmp_kmer(revcmp, it->first, k);
        fprint_kmer(outfile, revcmp, k);
        fprintf(outfile, "\t");
        if (right_ext == 'F' || right_ext == 'X')
            fprintf(outfile, "%c", right_ext);
        else
            fprintf(outfile, "%c", inv_base(right_ext));
        if (left_ext == 'F' || left_ext == 'X')
            fprintf(outfile, "%c", left_ext);
        else
            fprintf(outfile, "%c", inv_base(left_ext));
        fprintf(outfile, "\n");
    }
}

int main(int argc, char* argv[])
{
    mpi::environment env(argc, argv);
    mpi::communicator world;

    if (argc < 3) {
        printf("hammer outfile infile...\n");
        exit(1);
    }

    FastQReader* reader = get_reader(argc - 2, &argv[2], world, k);

    vector<string> prefix_boundaries;

    HashMap<kmer_t, qual_counts_t> kmer_map(INITIAL_CAPACITY, 0,
            (hash_map_hash_func_t) kmer_hash_K, (hash_map_eq_func_t) kmer_eq_K,
            kmer_size(k));
    build_map(reader, kmer_map, world);

    char outname[100];
    sprintf(outname, "%s.%d", argv[1], world.rank());
    FILE* outfile = fopen(outname, "w");
    distrib_print_ufx(outfile, kmer_map);
    fclose(outfile);
    //HashMap<kmer_t, ext_t> all_kmer_map(INITIAL_CAPACITY, 0,
    //        (hash_map_hash_func_t) kmer_hash_K, (hash_map_eq_func_t) kmer_eq_K,
    //        kmer_size(k));
    //gather_kmers(all_kmer_map, kmer_map, world);

    //if (world.rank() == 0) {
    //    FILE* outfile = fopen(argv[1], "w");
    //    print_ufxs(outfile, all_kmer_map);
    //    fclose(outfile);
    //}

    return 0;
}
