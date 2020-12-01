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
#include <sys/time.h>

// AVX Intrinsics (SIMD)
#include <immintrin.h>
// Threadpools: https://github.com/vit-vit/CTPL
#include "ctpl_stl.h"
#ifdef _OPENMP
#include <omp.h>
#endif

struct KrakenIndex {
  KrakenIndex(uint32_t c, uint32_t d) :
	count(c),
	dim(d),
	h_partitions(d, std::vector<float>()) { 
    this->dop = std::thread::hardware_concurrency();
    this->v_partitions = new std::vector<std::vector<float>>[this->dop];
  }

  KrakenIndex(uint32_t c, uint32_t d, uint32_t dp) :
	count(c),
	dim(d),
	h_partitions(d, std::vector<float>()),
	dop(dp) { 
    this->v_partitions = new std::vector<std::vector<float>>[this->dop];
  }

  ~KrakenIndex() {
  }

  uint32_t count;
  uint32_t dim;
  uint32_t dop;
  std::vector<std::vector<float> > h_partitions;
  std::vector<std::vector<float> >* v_partitions;
};

KrakenIndex* create_kraken(uint8_t dim);
KrakenIndex* create_kraken(uint8_t dim, uint32_t dop);
void insert(KrakenIndex* index, std::vector<float> point);
std::vector<uint32_t> range(KrakenIndex* index, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> range_simd(KrakenIndex* index, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> parallel_range(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> parallel_range_bitwise(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> parallel_range_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> parallel_range_bitmask(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> parallel_simd_scan(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint32_t> partitioned_range(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper, uint32_t dop);
std::vector<uint32_t> partitioned_range_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper, uint32_t dop);

#endif
