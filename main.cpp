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
#include "kmer_count_map.h"
#include "kmer_ext_map.h"
#include "contig.h"
#include "utils.h"
#include "lsh.h"
#include "contig_store.h"
#include "kmer_contig_map.h"

using namespace std;
namespace mpi = boost::mpi;
namespace fs  = boost::filesystem;


#define KMER_BUFFER_SIZE 1024
#define KMER_TAG 0
#define KMER_SIZE_TAG 1
#define DONE_TAG 2

#define LSH_ON 0

#ifndef LOAD_FROM_UFX
// Set to 0 if you want to run the kmer-count stage
// Set to 1 if you want to load from *.ufx.* instead of reading from fastq file
#define LOAD_FROM_UFX 1
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

void build_store(FastQReader* r, KmerCountMap& kmer_count_map, mpi::communicator& world)
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
                kmer_count_map.insert(recv_qekmer);
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

void print_ufxs(const char* outprefix, KmerExtMap& kmer_ext_map, int rank)
{
    stringstream ss;
    ss << outprefix << ".ufx." << rank;
    FILE* outfile = fopen(ss.str().c_str(), "w");
    if (outfile == NULL)
        panic("Could not open file: %s\n", ss.str().c_str());
    kmer_ext_map.print_ufxs(outfile);
    fclose(outfile);
}

void load_ufxs(char* file_prefix, KmerExtMap& kmer_ext_map, int rank)
{
    stringstream ss;
    ss << file_prefix << ".ufx." << rank;
    FILE* infile = fopen(ss.str().c_str(), "r");
    if (infile == NULL)
        panic("Could not open file: %s\n", ss.str().c_str());
    kmer_ext_map.load_ufxs(infile);
    fclose(infile);
}

void print_contigs(char* outprefix, ContigStore& contig_store, int rank)
{
    stringstream ss;
    ss << outprefix << ".contig." << rank;
    FILE* outfile = fopen(ss.str().c_str(), "w");
    if (outfile == NULL)
        panic("Could not open file: %s\n", ss.str().c_str());
    contig_store.fprint_contigs(outfile);
    fclose(outfile);
}

void print_contigs(char* outprefix, KmerContigMap& kmer_contig_map, int rank)
{
    stringstream ss;
    ss << outprefix << ".contig." << rank;
    FILE* outfile = fopen(ss.str().c_str(), "w");
    if (outfile == NULL)
        panic("Could not open file: %s\n", ss.str().c_str());
    kmer_contig_map.fprint_contigs(outfile);
    fclose(outfile);
}

/* Collects contigs on to one process (rank == 0). */
void gather_contigs(KmerContigMap& kmer_contig_map, ContigStore& contig_store, mpi::communicator& world)
{
    // TODO: Do we need to free all the kmers in the kmer_count_map?

    typedef struct {
        size_t size;
        base left_ext;
        base right_ext;
        char s[0];
    } __attribute__((packed)) contig_packet_t;

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
                contig_packet_t* cpacket;
                size_t size;
                while ((status = nethub.vrecv(i, (void**) &cpacket, &size)) == 0) {
                    Contig* contig = new Contig();
                    contig->left_ext = cpacket->left_ext;
                    contig->right_ext = cpacket->right_ext;
                    contig->s = std::string(cpacket->s, cpacket->size);
                    //contig->verify();
                    kmer_contig_map.insert(contig);
                    free(cpacket);
                }

                if (status == 1)
                    node_done[i] = true;
                else
                    all_done = false;

            }
        }
    }

    for (vector<Contig*>::iterator it = contig_store.contigs.begin();
            it != contig_store.contigs.end();
            it++) {
        Contig* c = *it;
        //c->verify();
        if (world.rank() == 0) {
            kmer_contig_map.insert(c);
        } else {
            size_t cpacket_size = sizeof(contig_packet_t) + c->s.size();
            contig_packet_t* cpacket = (contig_packet_t*) malloc(cpacket_size);
            cpacket->size = c->s.size();
            cpacket->left_ext = c->left_ext;
            cpacket->right_ext = c->right_ext;
            memcpy(cpacket->s, c->s.c_str(), cpacket->size);
            nethub.vsend(0, cpacket, cpacket_size);
            free(cpacket);
        }
    }
    nethub.done();
}

#if 0
void walk(Contig* contig, ContigStore& contig_store)
{
    kmer_a kmer[kmer_size(k)];
    bool revcmp_found = false;
    while(true) {
        contig->get_ext_kmer(kmer);
        ContigStore::iterator ext = contig_store.find(kmer, revcmp_found);
        if (contig_store.is_end(ext)) break;

        Contig* next_contig = ext->second;
        if (contig == next_contig) break;
        if (next_contig->s.size() == 0) break;

        if (revcmp_found) {
            next_contig->revcmp();
        }
        // FIXME: There are probably memory leak issues with not freeing kmers
        contig->s.append(next_contig->s, k-1, next_contig->s.size() - (k-1));
        contig->right_ext = next_contig->right_ext;

        next_contig->s.clear();
    }
}

void build_contigs(ContigStore& contig_store)
{
    for (ContigStore::iterator it = contig_store.begin();
            it != contig_store.end();
            it++)
    {
        Contig* contig = it->second;
        if (contig->s.size() == 0) continue;

        walk(contig, contig_store);
        contig->revcmp();
        walk(contig, contig_store);
        contig_store.add_to_final_contigs(it->second);
    }
    contig_store.trim();
}
#endif

int main(int argc, char* argv[])
{
    mpi::environment env(argc, argv);
    mpi::communicator world;

    if (argc < 3) {
        printf("hammer outfile infile...\n");
        exit(1);
    }

    //if (world.rank() == 0)
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
    KmerExtMap kmer_ext_map(k);

#if !LOAD_FROM_UFX
    {
        KmerCountMap kmer_count_map(k);
        FastQReader* reader = get_reader(argc - 2, &argv[2], world, k);
        build_store(reader, kmer_count_map, world);
        kmer_count_map.trim(kmer_ext_map);
    }
    print_ufxs(argv[1], kmer_ext_map, world.rank());
#else
    load_ufxs(argv[1], kmer_ext_map, world.rank());
#endif

    /* =======================
     * Phase 2: Contig walking
     * ======================= */
    ContigStore contig_store(k);
    KmerContigMap kmer_contig_map(k);
    kmer_ext_map.build_contigs(contig_store);

    //print_contigs(argv[1], contig_store, world.rank());

    gather_contigs(kmer_contig_map, contig_store, world);

    if (world.rank() == 0) {
        ContigStore joined_contig_store(k);
        kmer_contig_map.join_contigs(joined_contig_store);
        //build_contigs(contig_store);
        print_contigs(argv[1], joined_contig_store, world.rank() + 100);
    }

    return 0;
}
