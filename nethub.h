#ifndef _NETHUB_H_
#define _NETHUB_H_

#include <stdint.h>

#include <boost/mpi.hpp>

#include "kmer.h"

const size_t BUFFER_SIZE = 1024;

typedef struct {
    size_t size;
    uint8_t buffer[BUFFER_SIZE];
} nethub_packet_t;

class NetHub {
    enum packet_type {
        ACTIVE_TAG
    };

public:
    NetHub(boost::mpi::communicator& world, size_t data_size);
    ~NetHub();

    void send(size_t node_id, void* data);
    void done();

    /**
     * Returns:
     * -1 if there's nothing to read.
     * 0 if read was successful.
     * 1 if we got DONE_TAG.
     */
    int recv(size_t node_id, void* data);

private:
    void flush(size_t node_id);

    size_t data_size;

    boost::mpi::communicator& world;
    nethub_packet_t* send_packets;

    nethub_packet_t* recv_packets;
    size_t* recv_packet_idxs;
    boost::mpi::request* recv_requests;
};

#endif /* _NETHUB_H_ */
