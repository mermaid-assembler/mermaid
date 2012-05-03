#ifndef _HASH_MAP_H_
#define _HASH_MAP_H_

#include <cstdlib>

#include <sparsehash/sparse_hash_map>

#include "config.h"
#include "kmer.h"

using google::sparse_hash_map;

typedef size_t (*hash_map_hash_func_t)(size_t seed, void* key);
typedef bool (*hash_map_eq_func_t)(void* x, void* y);

/* NOTE - KeyType MUST be a pointer. */
template <class KeyType, class ValType>
class HashMap {
protected:
    class HashFunc {
    public:
        HashFunc(size_t seed, hash_map_hash_func_t hash) : seed(seed), hash(hash) { }
        size_t operator()(KeyType key) const { return hash(seed, key); }

    protected:
        size_t seed;
        hash_map_hash_func_t hash;
    };

    class EqFunc {
    public:
        EqFunc(hash_map_eq_func_t eq) : eq(eq) { }
        bool operator()(const KeyType& x, const KeyType& y) const { return eq(x, y); }

    protected:
        hash_map_eq_func_t eq;
    };

public:
    /* alloc_size is the size of the key for when new keys need to be added to
     * the map. */
    HashMap(size_t initial_capacity, size_t seed, hash_map_hash_func_t hash,
            hash_map_eq_func_t eq, size_t alloc_size)
        : initial_capacity(initial_capacity), seed(seed), hash_func(seed, hash),
        eq_func(eq), alloc_size(alloc_size),
         map(initial_capacity, hash_func, eq_func)
    { }

    // FIXME: Why does this silently fail? Can't we return NULL or a pointer
    // to the newly allocated ValType?
    /* Tries to insert key into the hash map. Fails silently if the key already
     * exists in the map. */
    void try_insert(KeyType key)
    {
        if (map.find(key) != map.end()) return;

        KeyType new_key = (KeyType) malloc(alloc_size);
        memcpy(new_key, key, alloc_size);
        map.insert(std::pair<const KeyType, ValType>(new_key, ValType()));
    }

    typedef sparse_hash_map<const KeyType, ValType, HashFunc, EqFunc> map_type_t;

protected:
    size_t initial_capacity;
    size_t seed;
    HashFunc hash_func;
    EqFunc eq_func;
    size_t alloc_size;

public:
    map_type_t map;
};

#endif /* _HASH_MAP_H_ */
