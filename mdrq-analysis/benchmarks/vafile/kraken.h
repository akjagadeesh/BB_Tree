#ifndef __Kraken_H
#define __Kraken_H

#include <cmath>

#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>

// VA file
#include "vafile.h"
// Threadpools: https://github.com/vit-vit/CTPL
#include "ctpl_stl.h"

struct KrakenIndex {
  KrakenIndex(uint32_t c, uint32_t d, uint32_t ddop) : count(c), dim(d), dop(ddop) {
//    this->dop = std::thread::hardware_concurrency();
	  this->partitions = new std::vector<std::vector<float>>[this->dop];
  }

  ~KrakenIndex() {
	  delete [] this->va_files;
	  delete [] this->partitions;
	  delete [] this->vafile_points;
  }

  uint32_t count;
  uint32_t dim;
  uint32_t dop;

  VAFile** va_files;
  std::vector<std::vector<float> >* vafile_points;
  std::vector<std::vector<float> >* partitions;
};

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop);
void partition_similar_rows(KrakenIndex* index);
void partition_similar_rows_different_partition(KrakenIndex* index);
void load_partitions(KrakenIndex* index);
void insert(KrakenIndex* index, std::vector<float> point);
void load_vafiles(KrakenIndex* index, std::vector<std::vector<float> > vafile_points);
std::vector<uint64_t> partitioned_range_vafile(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
std::vector<uint64_t> partitioned_range_vafile_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper);
#endif
