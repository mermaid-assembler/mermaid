#ifndef _TRIE_H_
#define _TRIE_H_

#include <cstdio>
#include <cstdlib>
#include <string>

#include "kmer.h"

class TrieNode {
public:
    TrieNode() { }
    void mark_ufx(FILE* ofile, kmer_t kmer, k_t depth, k_t max_depth) { }
};

class TrieBranchNode : public TrieNode {
public:
    TrieBranchNode();
    TrieNode* get_child(base b);
    void add_child(base b, TrieNode* child);
    void mark_ufx(FILE* ofile, kmer_t kmer, k_t depth, k_t max_depth);

protected:
    TrieNode* children[BASE::NUM_BASES];
};

class TrieLeafNode : public TrieNode {
public:
    TrieLeafNode();
    void mark_ufx(FILE* ofile, kmer_t kmer, k_t depth, k_t max_depth);

    count_t lcount[BASE::NUM_BASES];
    count_t rcount[BASE::NUM_BASES];
};

class Trie {
public:
    Trie(k_t k, std::string& outfile_name);
    /* TODO - Write a destructor which will free internal structures. */
    ~Trie();
    
    /* Here qkmer should be a k+2-mer. */
    void add_qkmer(qkmer_t* qkmer);

    void mark_ufx();

protected:
    k_t k;
    TrieBranchNode root;
    FILE* ofile;
};

#endif /* _TRIE_H_ */
