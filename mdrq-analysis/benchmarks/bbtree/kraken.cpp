#include "kraken.h"

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop) {
  return new KrakenIndex(0, dim, dop);
}

void insert(KrakenIndex* index, std::vector<float> point) {
  index->partitions[rand() % index->dop].push_back(point);
  index->count++;
}

void load_partitions(KrakenIndex* index) {
  index->kd_trees = new BBTree*[index->dop];
  for (size_t i = 0; i < index->dop; ++i)
    index->bb_trees[i] = new BBTree();

  for (size_t i = 0; i < index->dop; ++i) {
    for (size_t j = 0; j < index->partitions[i].size(); ++j) {
      index->bb_trees[i]->InsertObject(index->partitions[i][j], rand() % index->count);
    }
  }
}

std::vector<uint32_t> intersect(std::vector<uint32_t> first, std::vector<uint32_t> second) {
  std::vector<uint32_t> intersection;
  std::set_intersection(first.begin(), first.end(), second.begin(),
                        second.end(), std::back_inserter(intersection));

  return intersection;
}

inline void scan_partition_kdtree(int id, BBTree* bbtree, std::vector<uint32_t> &results, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> tmp_results = bbtree->SearchRange(lower, upper);
  results.resize(tmp_results.size());
  results = tmp_results;
}

inline void scan_partition_bbtree_simd(int id, BBTree* bbtree, std::vector<uint32_t> &results, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> tmp_results = bbtree->SearchRangeMT(lower, upper);
  results.resize(tmp_results.size());
  results = tmp_results;
}

std::vector<uint32_t> partitioned_range_bbtree(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(index->dop, std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[index->dop];

  for (uint32_t i = 0; i < index->dop; i++)
    futures[i] = pool->push(std::ref(scan_partition_bbtree), index->bb_trees[i], std::ref(intermediate_results[i]), lower, upper);

  for (uint32_t i = 0; i < index->dop; i++) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }
  delete [] futures;

  return results;
}

std::vector<uint32_t> partitioned_range_kdtree_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(index->dop, std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[index->dop];
  int partitions_visited = 0;

  for (uint32_t i = 0; i < index->dop; i++)
    futures[i] = pool->push(std::ref(scan_partition_kdtree_simd), index->bb_trees[i], std::ref(intermediate_results[i]), lower, upper);

  for (uint32_t i = 0; i < index->dop; i++) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }
  delete [] futures;

  return results;
}
