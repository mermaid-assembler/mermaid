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
#include "kmer_count_store.h"
#include "contig.h"
#include "utils.h"
#include "lsh.h"

using namespace std;
namespace mpi = boost::mpi;
namespace fs  = boost::filesystem;


#define KMER_BUFFER_SIZE 1024
#define KMER_TAG 0
#define KMER_SIZE_TAG 1
#define DONE_TAG 2

#define LSH_ON 0

#ifndef KMER_COUNT_ON
// Set to 0 if you want to load from *.ufx.* instead of reading from fastq file
// Set to 1 if you want to run the kmer-count stage
#define KMER_COUNT_ON 1
#endif

static const k_t k = K;

int get_kmer_bin(qekmer_t* qekmer, k_t k, int world_size)
{
#if LSH_ON
    size_t hash = lsh(qekmer->kmer, k);
#else
    size_t hash = kmer_hash(0, qekmer->kmer, k);
#endif
    return hash % world_size;
}

/* Checks to see if the qekmer meets the criteria to being sent. For example,
 * this function checks to make sure that the kmer does not contain an 'N'
 * anywhere in the kmer.
 */
bool check_qekmer_qual(qekmer_t* qekmer, k_t k)
{
    if (qekmer->lqual <= Q_MIN && qekmer->rqual <= Q_MIN)
        return false;

    return true;
}

/* Sets the qekmer to the canonical qekmer. */
void canonize_qekmer(qekmer_t* qekmer, k_t k)
{
    kmer_a revcmp[kmer_size(k)];
    revcmp_kmer(revcmp, qekmer->kmer, k);
    if (cmp_kmer(qekmer->kmer, revcmp, k) > 0) {
        memcpy(qekmer->kmer, revcmp, kmer_size(k));
        qual_t tmp_qual = qekmer->lqual;
        base tmp_ext = inv_base((base) qekmer->exts.left);
        qekmer->lqual = qekmer->rqual;
        qekmer->exts.left = inv_base((base) qekmer->exts.right);
        qekmer->rqual = tmp_qual;
        qekmer->exts.right = tmp_ext;
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

void build_store(FastQReader* r, KmerCountStore& kmer_store, mpi::communicator& world)
{
    NetHub nethub(world, qekmer_size(k));

    bool node_done[world.size()]; // Stores the non-blocking receive requests
    for (int i = 0; i < world.size(); i++) {
        node_done[i] = false;
    }
    bool done_reading = false;
    bool all_done = false;

    qekmer_t* send_qekmer = (qekmer_t*) malloc(qekmer_size(k));
    qekmer_t* recv_qekmer = (qekmer_t*) malloc(qekmer_size(k));

    while (!all_done) {
        if (r->read_next(send_qekmer)) {
            if (check_qekmer_qual(send_qekmer, k)) {
                canonize_qekmer(send_qekmer, k);
                int node_id = get_kmer_bin(send_qekmer, k, world.size());
                nethub.send(node_id, send_qekmer);
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
            while ((status = nethub.recv(i, recv_qekmer)) == 0) {
                kmer_store.insert(recv_qekmer);
            }

            if (status == 1)
                node_done[i] = true;
            else
                all_done = false;
        }
    }

    free(send_qekmer);
    free(recv_qekmer);
}

#if 0
void gather_kmers(HashMap<kmer_t, extensions_t>& all_kmer_map,
        HashMap<kmer_t, kmer_info_t>& kmer_map, mpi::communicator& world)
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

    HashMap<kmer_t, kmer_info_t>::map_type_t::iterator it = kmer_map.map.begin();
    while (!all_done) {
        if (it != kmer_map.map.end()) {
            const kmer_t& kmer = it->first;
            kmer_info_t& qual_counts = it->second;
            extensions_t extensions = qual_counts_2_extensions(&qual_counts);

            if (check_hq_extensions(extensions)) {
                if (world.rank() == 0) {
                    all_kmer_map.try_insert(kmer);
                    all_kmer_map.map[kmer] = extensions;
                } else {
                    for_base_in_kmer(b, kmer, k) {
                        set_base(send_ekmer->kmer, b_i_, b);
                    } end_for;
                    send_ekmer->ext = extensions;
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

void distrib_print_ufx(FILE* outfile, HashMap<kmer_t, kmer_info_t> kmer_map)
{
    for (HashMap<kmer_t, kmer_info_t>::map_type_t::iterator it = kmer_map.map.begin();
            it != kmer_map.map.end();
            it++) {
        const kmer_t& kmer = it->first;
        kmer_info_t& qual_counts = it->second;
        extensions_t extensions = qual_counts_2_extensions(&qual_counts);

        if (!check_hq_extensions(extensions))
            continue;

        char left_ext = 0;
        char right_ext = 0;

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (extensions.ext & 1 << (BASE::NUM_BASES + i)) {
                if (left_ext)
                    left_ext = 'F';
                else 
                    left_ext = base2char((base) i);
            }
        }
        if (!left_ext)
            left_ext = 'X';

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (extensions.ext & 1 << i) {
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

void print_ufxs(FILE* outfile, HashMap<kmer_t, extensions_t>& all_kmer_map)
{
    for (HashMap<kmer_t, extensions_t>::map_type_t::iterator it = all_kmer_map.map.begin();
            it != all_kmer_map.map.end();
            it++) {

        char left_ext = 0;
        char right_ext = 0;

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (it->second.ext & 1 << (BASE::NUM_BASES + i)) {
                if (left_ext)
                    left_ext = 'F';
                else 
                    left_ext = base2char((base) i);
            }
        }
        if (!left_ext)
            left_ext = 'X';

        for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
            if (it->second.ext & 1 << i) {
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

void build_contigs(vector<Contig*>& contigs, HashMap<kmer_t, kmer_info_t> kmer_map)
{
    for (HashMap<kmer_t, kmer_info_t>::map_type_t::iterator it = kmer_map.map.begin();
            it != kmer_map.map.end();
            it++) {
        const kmer_t& kmer = it->first;
        kmer_info_t& kmer_info = it->second;
        //extensions_t extensions = qual_counts_2_extensions(&kmer_info);

        //if (!check_hq_extensions(extensions))
        //    continue;
    }
}
#endif

void print_ufxs(const char* outprefix, KmerCountStore& kmer_store, int rank)
{
    stringstream ss;
    ss << outprefix << ".ufx." << rank;
    FILE* outfile = fopen(ss.str().c_str(), "w");
    if (outfile == NULL)
        panic("Could not open file: %s\n", ss.str().c_str());
    kmer_store.print_ufxs(outfile);
    fclose(outfile);
}

void load_ufxs(char* file_prefix, KmerCountStore& kmer_store, int rank)
{
    stringstream ss;
    ss << file_prefix << ".ufx." << rank;
    FILE* infile = fopen(ss.str().c_str(), "r");
    if (infile == NULL)
        panic("Could not open file: %s\n", ss.str().c_str());
    kmer_store.load_ufxs(infile);
    fclose(infile);
}

void print_contigs(char* outprefix, ContigStore& contig_store, int rank)
{
    stringstream ss;
    ss << outprefix << ".contig." << rank;
    FILE* outfile = fopen(ss.str().c_str(), "w");
    if (outfile == NULL)
        panic("Could not open file: %s\n", ss.str().c_str());
    contig_store.print_contigs(outfile);
    fclose(outfile);
}

/* Collects contigs on to one process (rank == 0). */
void gather_contigs(ContigStore& contig_store, mpi::communicator& world)
{
    // TODO: Do we need to free all the kmers in the kmer_count_store?

    typedef struct {
        size_t size;
        base left_ext;
        base right_ext;
        char s[0];
    } __attribute__((packed)) contig_info;

    NetHub nethub(world, 0);

    if (world.rank() == 0) {
        bool node_done[world.size()];
        for (int i = 0; i < world.size(); i++) {
            node_done[i] = false;
        }
        bool all_done = false;
        while (!all_done) {
            all_done = true;
            for (int i = 1; i < world.size(); i++) {
                if (node_done[i]) continue;

                int status;
                contig_info* cinfo;
                size_t size;
                while ((status = nethub.vrecv(i, (void**) &cinfo, &size)) == 0) {
                    Contig* contig = new Contig();
                    contig->left_ext = cinfo->left_ext;
                    contig->right_ext = cinfo->right_ext;
                    contig->s = std::string(cinfo->s, cinfo->size);
                    free(cinfo);
                    contig_store.add_contig(contig);
                }

                if (status == 1)
                    node_done[i] = true;
                else
                    all_done = false;

            }
        }
    }
    else /* world.rank() != 0 */
    {
        for (ContigStore::iterator it = contig_store.begin();
             it != contig_store.end();
             it++) {
            Contig* c = *it;
            size_t cinfo_size = sizeof(contig_info) + c->s.size();
            contig_info* cinfo = (contig_info*) malloc(cinfo_size);
            cinfo->size = c->s.size();
            cinfo->left_ext = c->left_ext;
            cinfo->right_ext = c->right_ext;
            memcpy(cinfo->s, c->s.c_str(), cinfo->size);
            nethub.vsend(0, cinfo, cinfo_size);
            free(cinfo);
        }
        nethub.done();
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

    //{
    //    printf("PID %d on %d ready for attach\n", getpid(), world.rank());
    //    fflush(stdout);
    //    volatile int i = 0;
    //    while (0 == i)
    //        sleep(5);
    //}

    /* =======================
     * Phase 1: k-mer counting
     * ======================= */
    KmerCountStore kmer_store(k);

#if KMER_COUNT_ON
    FastQReader* reader = get_reader(argc - 2, &argv[2], world, k);
    build_store(reader, kmer_store, world);
    kmer_store.trim();
    print_ufxs(argv[1], kmer_store, world.rank());
#else
    load_ufxs(argv[1], kmer_store, world.rank());
#endif

    /* =======================
     * Phase 2: Contig walking
     * ======================= */
    ContigStore contig_store;
    kmer_store.build_contigs(contig_store);

    print_contigs(argv[1], contig_store, world.rank());

    //gather_contigs(contig_store, world);

    //if (world.rank() == 0) {
    //    print_contigs(argv[1], contig_store, world.rank());
    //}

    //if (world.rank() == 0) {
    //    FILE* outfile = fopen(argv[1], "w");
    //    print_ufxs(outfile, all_kmer_map);
    //    fclose(outfile);
    //}

    return 0;
}
