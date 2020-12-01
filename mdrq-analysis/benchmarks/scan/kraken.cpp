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
#include "kraken.h"

static double gettime(void) {
  struct timeval now_tv;
  gettimeofday (&now_tv,NULL);
  return ((double)now_tv.tv_sec) + ((double)now_tv.tv_usec) / 1000000.0;
}

KrakenIndex* create_kraken(uint8_t dim) {
  return new KrakenIndex(0, dim);
}

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop) {
  return new KrakenIndex(0, dim, dop);
}

void insert(KrakenIndex* index, std::vector<float> point) {
  // insert into horizontal partitions
  for (uint8_t i = 0; i < index->dim; i++)
    index->h_partitions[i].push_back(point[i]);
  // insert into vertical partitions
  index->v_partitions[(rand() % index->dop)].push_back(point);
  index->count++;
}

bool compareTwoRows(float* a, float* b){
  return ( (a[0]<b[0]) || ((a[0]==b[0])&&(a[1]<b[1])) );
}

std::vector<uint32_t> intersect(std::vector<uint32_t> first, std::vector<uint32_t> second) {
  std::vector<uint32_t> intersection;
  std::set_intersection(first.begin(), first.end(), second.begin(),
                        second.end(), std::back_inserter(intersection));

  return intersection;
}

std::vector<uint32_t> intersect_iterators(std::vector<uint32_t>::iterator first_begin, std::vector<uint32_t>::iterator first_end,
                                          std::vector<uint32_t>::iterator second_begin, std::vector<uint32_t>::iterator second_end) {
  std::vector<uint32_t> intersection;
  std::set_intersection(first_begin, first_end, second_begin,
                        second_end, std::back_inserter(intersection));

  return intersection;
}

std::vector<uint32_t> range(KrakenIndex* index, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t>  results;
  for (size_t k = 0; k < index->dop; ++k) {
	  for (size_t i = 0; i < index->v_partitions[k].size(); ++i) {
	    bool match = true;
	    for (size_t j = 0; j < index->dim; ++j) {
	      if (index->v_partitions[k][i][j] < lower[j] || index->v_partitions[k][i][j] > upper[j]) {
		match = false; break;
	      }
	    }
	    if (match)
	      results.push_back(i);
	  }
  }

  return results;
}

std::vector<uint32_t> range_simd(KrakenIndex* index, std::vector<float> lower, std::vector<float> upper) {
  __m256 lower_reg, upper_reg, search_reg, lower_res, upper_res;
  size_t mask_lower, mask_upper, mask;

  size_t simd_compares = ((size_t) (index->dim / 8)) * 8;
  std::vector<uint32_t>  results;
  for (size_t k = 0; k < index->dop; ++k) {
	  for (size_t i = 0; i < index->v_partitions[k].size(); ++i) {
	    bool match = true;
	    size_t j = 0;

	    for (; j < simd_compares; j += 8) {
		    lower_reg = _mm256_loadu_ps(&lower[j]);
		    upper_reg = _mm256_loadu_ps(&upper[j]);
		    search_reg = _mm256_loadu_ps(&index->v_partitions[k][i][j]);
		    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
		    mask_lower = _mm256_movemask_ps(lower_res);
		    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
		    mask_upper = _mm256_movemask_ps(upper_res);
		    mask       = mask_lower & mask_upper;

		    // mask is set to 0xFF if every SIMD lane matchs
		    if (mask < 0xFF) {
			    match = false;
			    j = index->dim;
			    break;
		    }
	    }

	    for (; j < index->dim; ++j) {
		    if (lower[j] > index->v_partitions[k][i][j] ||
        	        upper[j] < index->v_partitions[k][i][j]) {
			    match = false; break;
		    }
	    }

	    if (match)
		    results.push_back(i);
	  }
  }

  return results;
}

inline void scan_dimension(int id, KrakenIndex* index, std::vector<uint32_t> &results, uint8_t dim, float lower, float upper) {
  for (size_t i = 0; i < index->count; ++i)
    if (index->h_partitions[dim][i] >= lower && index->h_partitions[dim][i] <= upper)
      results.push_back(i);
}

inline void scan_partition(int id, KrakenIndex* index, std::vector<uint32_t> &results, std::vector<float> lower, std::vector<float> upper, uint32_t partition) {
  for (size_t i = 0; i < index->v_partitions[partition].size(); ++i) {
    bool match = true;
    for (size_t j = 0; j < index->dim; ++j) {
      if (index->v_partitions[partition][i][j] < lower[j] || index->v_partitions[partition][i][j] > upper[j]) {
        match = false; break;
      }
    }
    if (match)
      results.push_back(i);
  }
}

inline void scan_partition_simd(int id, KrakenIndex* index, std::vector<uint32_t> &results, std::vector<float> lower, std::vector<float> upper, uint32_t partition) {
  __m256 lower_reg, upper_reg, search_reg, lower_res, upper_res;
  size_t mask_lower, mask_upper, mask;

  size_t simd_compares = ((size_t) (index->dim / 8)) * 8;
  for (size_t i = 0; i < index->v_partitions[partition].size(); ++i) {
    bool match = true;
    size_t j = 0;

    for (; j < simd_compares; j += 8) {
      lower_reg = _mm256_loadu_ps(&lower[j]);
      upper_reg = _mm256_loadu_ps(&upper[j]);
      search_reg = _mm256_loadu_ps(&index->v_partitions[partition][i][j]);
      lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
      mask_lower = _mm256_movemask_ps(lower_res);
      upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
      mask_upper = _mm256_movemask_ps(upper_res);
      mask       = mask_lower & mask_upper;

      // mask is set to 0xFF if every SIMD lane matchs
      if (mask < 0xFF) {
        match = false;
        j = index->dim;
        break;
      }
    }

    for (; j < index->dim; ++j) {
      if (lower[j] > index->v_partitions[partition][i][j] ||
	  upper[j] < index->v_partitions[partition][i][j]) {
	match = false; break;
      }
    }

    if (match)
      results.push_back(i);
  }
}

inline void scan_dimension_simd(int id, KrakenIndex* index, std::vector<uint32_t> &results, uint8_t dim, float lower, float upper) {
  if (lower == std::numeric_limits<float>::min() && upper == std::numeric_limits<float>::max())
    return;

  __m256 lower_reg = _mm256_set1_ps(lower);
  __m256 upper_reg = _mm256_set1_ps(upper);
  __m256 search_reg, lower_res, upper_res;
  size_t mask_lower, mask_upper, mask;

  uint32_t curPos = 0;
  uint32_t limit  = index->count;
  while (curPos < limit) {
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;

    if (mask & 0x1)
      results.push_back(curPos);
    if (mask & 0x2)
      results.push_back(curPos + 1);
    if (mask & 0x4)
      results.push_back(curPos + 2);
    if (mask & 0x8)
      results.push_back(curPos + 3);
    if (mask & 0x10)
      results.push_back(curPos + 4);
    if (mask & 0x20)
      results.push_back(curPos + 5);
    if (mask & 0x40)
      results.push_back(curPos + 6);
    if (mask & 0x80)
      results.push_back(curPos + 7);
    curPos += 8;
  }
}

inline void scan_dimension_simd_bitmask(int id, KrakenIndex* index, std::vector<uint64_t> &results, uint8_t dim, float lower, float upper) {
  __m256 lower_reg = _mm256_set1_ps(lower);
  __m256 upper_reg = _mm256_set1_ps(upper);
  __m256 search_reg, lower_res, upper_res;
  uint64_t mask_lower, mask_upper, mask;

  uint32_t curPos = 0;
  uint32_t limit  = index->count;
  while (curPos < limit) {
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    // mask is > 0x00 if a match has been found
    if (mask)
      results[curPos / 64] = (mask << 56);
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos + 8]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    if (mask)
      results[curPos / 64] |= (mask << 48);
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos + 16]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    if (mask)
      results[curPos / 64] |= (mask << 40);
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos + 24]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    if (mask)
      results[curPos / 64] |= (mask << 32);
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos + 32]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    if (mask)
      results[curPos / 64] |= (mask << 24);
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos + 40]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    if (mask)
      results[curPos / 64] |= (mask << 16);
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos + 48]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    if (mask)
      results[curPos / 64] |= (mask << 8);
    search_reg = _mm256_loadu_ps(&index->h_partitions[dim][curPos + 56]);
    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
    mask_lower = _mm256_movemask_ps(lower_res);
    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
    mask_upper = _mm256_movemask_ps(upper_res);
    mask       = mask_lower & mask_upper;
    if (mask)
      results[curPos / 64] |= mask;
    curPos += 64;
  }
}

inline void scan_dimension_bitmask(int id, KrakenIndex* index, std::vector<uint64_t> &results, uint8_t dim, float lower, float upper) {
  uint32_t bitmasks = (index->count / 64);

  for (size_t i = 0; i < bitmasks; ++i) {
    for (size_t j = 0; j < 64; ++j) {
      results[i] |= (((uint64_t) (index->h_partitions[dim][i * 64 + j] >= lower &&
                      index->h_partitions[dim][i * 64 + j] <= upper))
                     << 63 - j);
    }
  }
}

std::vector<uint32_t> partitioned_range(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper, uint32_t dop) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(dop, std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[dop];

  for (uint32_t i = 0; i < dop; ++i) {
    futures[i] = pool->push(std::ref(scan_partition), index, std::ref(intermediate_results[i]), lower, upper, i);
  }

  for (uint32_t i = 0; i < dop; ++i) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }
  delete [] futures;

  return results;
}

std::vector<uint32_t> partitioned_range_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper, uint32_t dop) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(dop, std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[dop];

  uint32_t step = std::ceil(index->count/(float)dop);
  for (uint32_t i = 0; i < dop; i++)
    futures[i] = pool->push(std::ref(scan_partition_simd), index, std::ref(intermediate_results[i]), lower, upper, i);

  for (uint32_t i = 0; i < dop; i++) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }
  delete [] futures;

  return results;
}

std::vector<uint32_t> parallel_range(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  double start = gettime();
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(
                             index->dim,std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[index->dim];

  for (size_t i = 0; i < index->dim; ++i) {
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max()) {
      continue;
    }
    futures[i] = pool->push(std::ref(scan_dimension), index,
           		    std::ref(intermediate_results[i]),
	         	    i, lower[i], upper[i]);
  }

  size_t dim = 0;
  for (size_t i = 0; i < index->dim; ++i) {
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max()) {
      continue;
    }
    futures[i].get();
    results = intermediate_results[i];
    dim = i + 1;
    break;
  }

  for (size_t i = dim; i < index->dim; ++i) {
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    futures[i].get();
    results = intersect(results, intermediate_results[i]);
  }
  delete [] futures;

  return results;
}

std::vector<uint32_t> parallel_range_bitwise(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  uint32_t bitmask_count = index->count/64;
  std::vector<uint32_t> tid_results;
  std::vector<uint64_t> results(bitmask_count, 0x0000000000000000);
  std::vector<std::vector<uint64_t>> bitmasks(index->dim, std::vector<uint64_t>(bitmask_count, 0x0000000000000000));
  std::future<void> *futures = new std::future<void>[index->dim];

  for (uint8_t i = 0; i < index->dim; ++i) {
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    futures[i] = pool->push(std::ref(scan_dimension_bitmask), index, std::ref(bitmasks[i]), i, lower[i], upper[i]);
  }

  size_t dim = 0;
  for (size_t i = 0; i < index->dim; ++i) {
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    futures[i].get();
    results = bitmasks[i];
    dim = i + 1;
    break;
  }

  for (size_t j = dim; j < index->dim; ++j) {
    if (lower[j] == std::numeric_limits<float>::min() &&
        upper[j] == std::numeric_limits<float>::max())
      continue;
    futures[j].get();
    for(uint32_t i = 0; i < bitmask_count; ++i)
      results[i] &= bitmasks[j][i];
  }
  delete [] futures;

  for (size_t i = 0; i < bitmask_count; ++i)
    if (results[i] > 0x0000000000000000)
      for (size_t j = 0; j < 64; ++j)
        if ((results[i] >> j) & 0x1)
	  tid_results.push_back((uint32_t) (i*64 + (63-j)));

  return tid_results;
}

std::vector<uint32_t> parallel_range_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(index->dim,std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[index->dim];

  for (uint8_t i = 0; i < index->dim; i++)
    futures[i] = pool->push(std::ref(scan_dimension_simd), index, std::ref(intermediate_results[i]), i, lower[i], upper[i]);

  size_t dim = 0;
  for (size_t i = 0; i < index->dim; ++i) {
    futures[i].get();
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    results = intermediate_results[i];
    dim = i + 1;
    break;
  }

  for (size_t i = dim; i < index->dim; ++i) {
    futures[i].get();
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    results = intersect(results, intermediate_results[i]);
  }
  delete [] futures;

  return results;
}

std::vector<uint32_t> parallel_range_bitmask(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(index->dim,std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[index->dim];

  for (uint8_t i = 0; i < index->dim; i++) {
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    futures[i] = pool->push(std::ref(scan_dimension_simd), index, std::ref(intermediate_results[i]), i, lower[i], upper[i]);
  }

  size_t dim = 0;
  for (size_t i = 0; i < index->dim; ++i) {
    futures[i].get();
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    results = intermediate_results[i];
    dim = i + 1;
    break;
  }

  for (size_t i = dim; i < index->dim; ++i) {
    futures[i].get();
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    results = intersect(results, intermediate_results[i]);
  }
  delete [] futures;

  return results;
}

std::vector<uint32_t> parallel_simd_scan(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  uint32_t bitmask_count = index->count/64;
  std::vector<uint32_t> tid_results;
  std::vector<uint64_t> results(bitmask_count, 0x0000000000000000);
  std::vector<std::vector<uint64_t>> bitmasks(index->dim, std::vector<uint64_t>(bitmask_count, 0x0000000000000000));
  std::future<void> *futures = new std::future<void>[index->dim];

  for (uint8_t i = 0; i < index->dim; ++i)
    futures[i] = pool->push(std::ref(scan_dimension_simd_bitmask), index, std::ref(bitmasks[i]), i, lower[i], upper[i]);

  size_t dim = 0;
  for (size_t i = 0; i < index->dim; ++i) {
    futures[i].get();
    if (lower[i] == std::numeric_limits<float>::min() &&
        upper[i] == std::numeric_limits<float>::max())
      continue;
    results = bitmasks[i];
    dim = i + 1;
    break;
  }

  for (size_t j = dim; j < index->dim; ++j) {
    futures[j].get();
    if (lower[j] == std::numeric_limits<float>::min() &&
        upper[j] == std::numeric_limits<float>::max())
      continue;
    for(uint32_t i = 0; i < bitmask_count; ++i)
      results[i] &= bitmasks[j][i];
  }
  delete [] futures;

  for (size_t i = 0; i < bitmask_count; ++i)
    if (results[i] > 0x0000000000000000)
      for (size_t j = 0; j < 64; ++j)
        if ((results[i] >> j) & 0x1)
	  tid_results.push_back((uint32_t) (i*64 + (63-j)));

  return tid_results;
}
