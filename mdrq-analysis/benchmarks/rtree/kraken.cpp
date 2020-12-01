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

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop) {
  return new KrakenIndex(0, dim, dop);
}

void insertIntoRTree(Index* rtree, std::vector<float> point, uint32_t id) {
  double* coords = new double[point.size()];
  for (size_t i = 0; i < point.size(); i++)
    coords[i] = (double) point[i];
  uint8_t* pData = 0;
  uint32_t nDataLength = 0;

  SpatialIndex::IShape* shape = new SpatialIndex::Point(coords, point.size());
  rtree->index().insertData(nDataLength,pData,*shape,id+1);
  delete shape;
}

void insert(KrakenIndex* index, std::vector<float> point) {
  index->partitions[rand() % index->dop].push_back(point);
  //insertIntoRTree(index->rtrees[rand() % index->dop], point, index->count);
  index->count++;
}

void load_rtrees(KrakenIndex* index) {
  for (size_t i = 0; i < index->dop; ++i) {
    for (size_t j = 0; j < index->partitions[i].size(); ++j) {
      insertIntoRTree(index->rtrees[i],
		      index->partitions[i][j],
		      rand() % index->count);
    }
  }
}

inline void scan_partition_rtree(int id, Index* rtree, std::vector<uint32_t> &results, std::vector<float> lower, std::vector<float> upper) {
  uint32_t dim = lower.size();
  double* plow  = new double[dim];
  double* phigh = new double[dim];
  for (size_t i = 0; i < dim; i++) {
    plow[i]  = lower[i];
    phigh[i] = upper[i];
  }

  // get nearest maxResults shapes form index
  ObjVisitor* visitor = new ObjVisitor;
  SpatialIndex::Region* r = new SpatialIndex::Region(plow, phigh, dim);
  rtree->index().intersectsWithQuery(*r, *visitor);

  int64_t nResultCount = visitor->GetResultCount();

  // get actual results
  std::vector<SpatialIndex::IData*>& rresults = visitor->GetResults();

  // copy the Items into the newly allocated vector array
  // we need to make sure to clone the actual Item instead
  // of just the pointers, as the visitor will nuke them
  // upon destroy
  for (int64_t i = 0; i < nResultCount; i++)
    results.push_back(rresults[i]->getIdentifier());

  delete [] plow;
  delete [] phigh;
  delete r;
  delete visitor;
}

inline void scan_partition_rtree_simd(int id, Index* rtree, std::vector<uint32_t> &results, std::vector<float> lower, std::vector<float> upper) {
  uint32_t dim = lower.size();
  double* plow  = new double[dim];
  double* phigh = new double[dim];
  for (size_t i = 0; i < dim; i++) {
    plow[i]  = lower[i];
    phigh[i] = upper[i];
  }

  // get nearest maxResults shapes form index
  ObjVisitor* visitor = new ObjVisitor;
  SpatialIndex::Region* r = new SpatialIndex::Region(plow, phigh, dim);
  rtree->index().intersectsWithQuerySIMD(*r, *visitor);

  int64_t nResultCount = visitor->GetResultCount();

  // get actual results
  std::vector<SpatialIndex::IData*>& rresults = visitor->GetResults();

  // copy the Items into the newly allocated vector array
  // we need to make sure to clone the actual Item instead
  // of just the pointers, as the visitor will nuke them
  // upon destroy
  for (int64_t i = 0; i < nResultCount; i++)
    results.push_back(rresults[i]->getIdentifier());

  delete [] plow;
  delete [] phigh;
  delete r;
  delete visitor;
}

// V_SCAN with R*-tree
std::vector<uint32_t> partitioned_range_rtree(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(index->dop, std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[index->dop];

  for (uint32_t i = 0; i < index->dop; i++)
    futures[i] = pool->push(std::ref(scan_partition_rtree), index->rtrees[i], std::ref(intermediate_results[i]), lower, upper);

  for (uint32_t i = 0; i < index->dop; i++) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }
  std::sort(std::begin(results), std::end(results));

  delete [] futures;

  return results;
}

// V_SCAN with R*-tree and SIMD
std::vector<uint32_t> partitioned_range_rtree_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint32_t> results;
  std::vector<std::vector<uint32_t>> intermediate_results(index->dop, std::vector<uint32_t>());
  std::future<void> *futures = new std::future<void>[index->dop];

  for (uint32_t i = 0; i < index->dop; i++)
    futures[i] = pool->push(std::ref(scan_partition_rtree_simd), index->rtrees[i], std::ref(intermediate_results[i]), lower, upper);

  for (uint32_t i = 0; i < index->dop; i++) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }
  std::sort(std::begin(results), std::end(results));

  delete [] futures;

  return results;
}
