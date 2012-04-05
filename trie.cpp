#include <cassert>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <string>

#include "trie.h"
#include "kmer.h"
#include "config.h"

TrieBranchNode::TrieBranchNode()
    : TrieNode()
{
    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        children[i] = NULL;
    }
}

TrieNode* TrieBranchNode::get_child(base b)
{
    return children[b];
}

void TrieBranchNode::add_child(base b, TrieNode* child)
{
    children[b] = child;
}

// kmer is modified as we traverse. Must traverse in a DFS order.
void TrieBranchNode::mark_ufx(FILE* ofile, kmer_t kmer, k_t depth, k_t max_depth)
{
    //assert(depth < max_depth);
    for (uint8_t i = 0; i < BASE::NUM_BASES; i++)
    {
        if (children[i] != NULL) {
            set_base(kmer, depth, i);
            if (depth == max_depth - 1) {
                static_cast<TrieLeafNode*>(children[i])->mark_ufx(ofile, kmer, depth + 1, max_depth);
            } else {
                static_cast<TrieBranchNode*>(children[i])->mark_ufx(ofile, kmer, depth + 1, max_depth);
            }
        }
    }
}

TrieLeafNode::TrieLeafNode()
    : TrieNode()
{
    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        lcount[i] = 0;
        rcount[i] = 0;
    }
}

void TrieLeafNode::mark_ufx(FILE* ofile, kmer_t kmer, k_t depth, k_t max_depth)
{
    assert(depth == max_depth);

    char left_ext = 0;
    char right_ext = 0;

    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        if (lcount[i] >= D_MIN) {
            if (left_ext)
                left_ext = 'F';
            else
                left_ext = base2char((base) i);
        }
    }

    if (!left_ext)
        left_ext = 'X';

    for (uint8_t i = 0; i < BASE::NUM_BASES; i++) {
        if (rcount[i] >= D_MIN) {
            if (right_ext)
                right_ext = 'F';
            else
                right_ext = base2char((base) i);
        }
    }

    if (!right_ext)
        right_ext = 'X';

    if (left_ext != 'X' || right_ext != 'X') {
        fprint_kmer(ofile, kmer, depth);
        fprintf(ofile, "\t%c%c\n", left_ext, right_ext);

        uint8_t revcmp[kmer_size(depth)];
        revcmp_kmer(revcmp, kmer, depth);
        fprint_kmer(ofile, revcmp, depth);
        fprintf(ofile, "\t");
        if (right_ext == 'F' || right_ext == 'X')
            fprintf(ofile, "%c", right_ext);
        else
            fprintf(ofile, "%c", inv_base(right_ext));
        if (left_ext == 'F' || left_ext == 'X')
            fprintf(ofile, "%c", left_ext);
        else
            fprintf(ofile, "%c", inv_base(left_ext));
        fprintf(ofile, "\n");
    }
}

Trie::Trie(k_t k, std::string& outfile_name)
    : k(k), root(), ofile(fopen(outfile_name.c_str(), "w"))
{
}

Trie::~Trie()
{
    // TODO clean up more stuff
    fclose(ofile);
}


void Trie::add_qkmer(qkmer_t* qkmer)
{
    TrieNode* current = &root;
    k_t i;
    for (i = 1; i < k+1; i++) {
        base b = get_base(qkmer->kmer, i);
        TrieNode* child = static_cast<TrieBranchNode*>(current)->get_child(b);
        if (child == NULL) {
            break;
        } else {
            current = child;
        }
    }

    if (i < k+1) {
        for (; i < k+1; i++) {
            base b = get_base(qkmer->kmer, i);
            TrieNode* child;
            if (i == k) {
                child = new TrieLeafNode();
            } else {
                child = new TrieBranchNode();
            }
            static_cast<TrieBranchNode*>(current)->add_child(b, child);
            current = child;
        }
    }

    if (qkmer->lqual > Q_MIN)
        static_cast<TrieLeafNode*>(current)->lcount[get_base(qkmer->kmer, 0)]++;
    if (qkmer->rqual > Q_MIN)
        static_cast<TrieLeafNode*>(current)->rcount[get_base(qkmer->kmer, k + 1)]++;
}

void Trie::mark_ufx()
{
    uint8_t kmer[kmer_size(k)];
    root.mark_ufx(ofile, kmer, 0, k);
}
