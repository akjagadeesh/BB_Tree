/*********************************************************
*
*  Research Work of Stefan Sprenger
*  https://www2.informatik.hu-berlin.de/~sprengsz/
*
*  Used solely for scholastic work in course CSCE 614 for
*  the course research project. Adaptations and additions
*  are marked with //ADDED ... //ADDED.
*  
*********************************************************/
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
#include "kdtree/Tree.h"
// Threadpools: https://github.com/vit-vit/CTPL
#include "ctpl_stl.h"

struct KrakenIndex {
  KrakenIndex(uint32_t c, uint32_t d, uint32_t dopp) : count(c), dim(d), dop(dopp) {
//    this->dop = std::thread::hardware_concurrency();
    this->partitions = new std::vector<std::vector<float>>[this->dop];
  }

  ~KrakenIndex() {
    delete [] kd_trees;
    delete [] this->partitions;
  }

  uint32_t count;
  uint32_t dim;
  uint32_t dop;

  // kd tree
  Tree** kd_trees;
  std::vector<std::vector<float> >* kdtree_points;
  std::vector<std::vector<float> >* partitions;
};

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop);
void insert(KrakenIndex* index, std::vector<float> point);
void load_partitions(KrakenIndex* index);
void load_kdtrees(KrakenIndex* index, std::vector<std::vector<float> > kdtree_points);
std::vector<uint32_t> partitioned_range_kdtree(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> partitioned_range_kdtree_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
#endif
