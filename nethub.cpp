#include <boost/mpi.hpp>

#include "nethub.h"

namespace mpi = boost::mpi;

NetHub::NetHub(mpi::communicator& world, size_t data_size)
    : data_size(data_size), world(world), send_packets(), recv_packets(),
    recv_packet_idxs(), recv_requests()
{
    send_packets = (nethub_packet_t*) malloc(world.size() * sizeof(nethub_packet_t));

    recv_packets = (nethub_packet_t*) malloc(world.size() * sizeof(nethub_packet_t));
    recv_packet_idxs = (size_t*) malloc(world.size() * sizeof(size_t));
    //recv_requests = (mpi::request*) malloc(world.size() * sizeof(mpi::request));
    recv_requests = new mpi::request[world.size()];

    for (int i = 0; i < world.size(); i++) {
        send_packets[i].size = 0;

        recv_packets[i].size = 0;
        recv_packet_idxs[i] = 0;
        recv_requests[i] =  world.irecv(i, ACTIVE_TAG, 
                (uint8_t*) &recv_packets[i], sizeof(nethub_packet_t));
    }
}

NetHub::~NetHub()
{
    delete[] recv_requests;
    free(recv_packet_idxs);
    free(recv_packets);
    free(send_packets);
}

void NetHub::send(size_t node_id, void* data)
{
    size_t size = send_packets[node_id].size;
    memcpy(&send_packets[node_id].buffer[size], data, data_size);
    send_packets[node_id].size += data_size;

    if ((send_packets[node_id].size + data_size) > BUFFER_SIZE) {
        flush(node_id);
    }
}

void NetHub::flush(size_t node_id)
{
    if (!send_packets[node_id].size) {
        return;
    }

    world.send(node_id, ACTIVE_TAG, (uint8_t*) &send_packets[node_id],
            sizeof(nethub_packet_t));
    send_packets[node_id].size = 0;
}


void NetHub::done()
{
    for (int i = 0; i < world.size(); i++) {
        flush(i);
        send_packets[i].size = 0;
        world.send(i, ACTIVE_TAG, (uint8_t*) &send_packets[i],
                sizeof(nethub_packet_t));
    }
}

int NetHub::recv(size_t node_id, void* data)
{
    /* Check to see if we have any request objects which have received packets.
     * If not, return -1 right away (non-blocking). */
    if (!recv_requests[node_id].test()) {
        return -1;
    }

    /* If we received a packet and the size is 0, that means the node that sent
     * that packet is done sending packets. */
    if (!recv_packets[node_id].size) {
        return 1;
    }

    size_t idx = recv_packet_idxs[node_id];
    memcpy(data, &recv_packets[node_id].buffer[idx], data_size);
    recv_packet_idxs[node_id] += data_size;

    if (recv_packet_idxs[node_id] >= recv_packets[node_id].size) {
        recv_requests[node_id] = world.irecv(node_id, ACTIVE_TAG,
                (uint8_t*) &recv_packets[node_id], sizeof(nethub_packet_t));
        recv_packet_idxs[node_id] = 0;
    }
    
    return 0;
}
