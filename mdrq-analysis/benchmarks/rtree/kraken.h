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
// R*tree
#include <spatialindex/capi/sidx_api.h>
#include <spatialindex/capi/sidx_impl.h>
#include <spatialindex/capi/sidx_config.h>
// Threadpools: https://github.com/vit-vit/CTPL
#include "ctpl_stl.h"

#ifdef _OPENMP
#include <omp.h>
#endif

struct KrakenIndex {
  KrakenIndex(uint32_t c, uint32_t d, uint32_t dop) : count(c), dim(d), dop(dop) {
    // initialize R*-trees
    this->rtrees = new Index*[this->dop];
    Tools::PropertySet* ps = GetDefaults();
    Tools::Variant var;
    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = RT_RTree;
    ps->setProperty("IndexType", var);
    var.m_val.ulVal = dim;
    ps->setProperty("Dimension", var);
    var.m_val.ulVal = RT_Memory;
    ps->setProperty("IndexStorageType", var);
    var.m_val.ulVal = 96;
    ps->setProperty("IndexCapacity", var);
    var.m_val.ulVal = 96;
    ps->setProperty("LeafCapacity", var);
    for (uint8_t i = 0; i < this->dop; i++)
      this->rtrees[i] = new Index(*ps);
    delete ps;

    this->partitions = new std::vector<std::vector<float>>[this->dop];
  }

  ~KrakenIndex() {
    delete [] rtrees;
    delete [] this->partitions;
  }

  uint32_t count;
  uint32_t dim;
  uint32_t dop;

  std::vector<std::vector<float> >* partitions;

  // R*tree
  Index** rtrees;
};

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop);
void load_rtrees(KrakenIndex* index);
void insert(KrakenIndex* index, std::vector<float> point);
std::vector<uint32_t> partitioned_range_rtree(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> partitioned_range_rtree_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);

#endif
