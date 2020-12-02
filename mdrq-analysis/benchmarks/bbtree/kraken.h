#ifndef __Kraken_H
#define __Kraken_H

#include <cmath>

#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>

// AVX Intrinsics (SIMD)
#include <immintrin.h>
// kd tree
#include "bbtree/BBTree.h"
// Threadpools: https://github.com/vit-vit/CTPL
#include "ctpl_stl.h"

struct KrakenIndex {
  KrakenIndex(uint32_t c, uint32_t d, uint32_t dopp) : count(c), dim(d), dop(dopp) {
//    this->dop = std::thread::hardware_concurrency();
    this->partitions = new std::vector<std::vector<float>>[this->dop];
  }

  ~KrakenIndex() {
    delete [] bb_trees;
    delete [] this->partitions;
  }

  uint32_t count;
  uint32_t dim;
  uint32_t dop;

  // bb tree
  BBTree** bb_trees;
  std::vector<std::vector<float> >* bbtree_points;
  std::vector<std::vector<float> >* partitions;
};

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop);
void insert(KrakenIndex* index, std::vector<float> point);
void load_partitions(KrakenIndex* index);
void load_bbtrees(KrakenIndex* index, std::vector<std::vector<float> > bbtree_points);
std::vector<uint32_t> partitioned_range_bbtree(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> partitioned_range_bbtree_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
#endif
