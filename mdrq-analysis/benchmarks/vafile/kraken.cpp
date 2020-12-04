#include "kraken.h"

KrakenIndex* create_kraken(uint8_t dim, uint32_t dop) {
  return new KrakenIndex(0, dim, dop);
}

void insert(KrakenIndex* index, std::vector<float> point) {
  index->partitions[rand() % index->dop].push_back(point);
}

float euclidean_distance(std::vector<float> first, std::vector<float> second) {
  float distance = 0.0;
  for (size_t i = 0; i < first.size(); ++i)
    distance += (second[i] - first[i]) * (second[i] - first[i]);
  distance = sqrt(distance);

  return distance;
}

void partition_similar_rows(KrakenIndex* index) {
  bool partitioning_changed = true;
  while (partitioning_changed) {
    partitioning_changed = false;
    std::vector<std::vector<float>> mean_values = std::vector<std::vector<float>>(index->dop, std::vector<float>(index->dim, 0.0));

    #pragma omp parallel for
    for (size_t partition = 0; partition < index->dop; ++partition) {
      for (size_t i = 0; i < index->partitions[partition].size(); ++i)
        for (size_t j = 0; j < index->dim; ++j)
          mean_values[partition][j] += index->partitions[partition][i][j];
      for (size_t j = 0; j < index->dim; ++j)
        mean_values[partition][j] = mean_values[partition][j] / index->partitions[partition].size();
    }

    for (size_t partition = 0; partition < index->dop; ++partition) {
      for (size_t i = 0; i < index->partitions[partition].size(); ++i) {
        bool   local_changed     = false;
        size_t nearest_partition = partition;
        float  distance          = euclidean_distance(mean_values[partition], index->partitions[partition][i]);
        for (size_t j = 0; j < index->dop; ++j) {
          if (euclidean_distance(mean_values[j], index->partitions[partition][i]) < distance) {
            distance = euclidean_distance(mean_values[j], index->partitions[partition][i]);
            nearest_partition = j;
            local_changed = partitioning_changed = true;
          }
        }
        if (local_changed) {
            index->partitions[nearest_partition].push_back(index->partitions[partition][i]);
            index->partitions[partition].erase(index->partitions[partition].begin() + i);
        }
      }
    }
  }
}

void partition_similar_rows_different_partition(KrakenIndex* index) {
  bool partitioning_changed = true;
  size_t num_clusters = (index->count / index->dop) + 1;
  size_t cluster_size = index->dop;
  size_t iterations = 0;
  std::vector<std::vector<std::vector<float>>> clusters =
          std::vector<std::vector<std::vector<float>>>(num_clusters);

  // init clusters
  size_t cluster = 0;
//  #pragma omp parallel for
  for (size_t i = 0; i < index->dop; ++i) {
          for (size_t j = 0; j < index->partitions[i].size(); ++j) {
                  clusters[cluster].push_back(index->partitions[i][j]);
                  if (clusters[cluster].size() == cluster_size)
                          cluster++;
          }
  }
std::cout << "Initialized clusters." << std::endl;

  while (partitioning_changed && iterations++ < 20) {
    partitioning_changed = false;
    std::vector<std::vector<float>> mean_values = std::vector<std::vector<float>>(num_clusters, std::vector<float>(index->dim, 0.0));

    #pragma omp parallel for
    for (cluster = 0; cluster < num_clusters; ++cluster) {
      for (size_t i = 0; i < clusters[cluster].size(); ++i)
        for (size_t j = 0; j < index->dim; ++j)
          mean_values[cluster][j] += clusters[cluster][i][j];
      for (size_t j = 0; j < index->dim; ++j)
        mean_values[cluster][j] = mean_values[cluster][j] / clusters[cluster].size();
    }
std::cout << "Determined mean values." << std::endl;

    for (cluster = 0; cluster < num_clusters; ++cluster) {
      for (size_t i = 0; i < clusters[cluster].size(); ++i) {
        bool   local_changed     = false;
        size_t nearest_partition = cluster;
        float  distance          = euclidean_distance(mean_values[cluster], clusters[cluster][i]);
        for (size_t j = 0; j < num_clusters; ++j) {
          if (euclidean_distance(mean_values[j], clusters[cluster][i]) < distance) {
            distance = euclidean_distance(mean_values[j], clusters[cluster][i]);
            nearest_partition = j;
            local_changed = partitioning_changed = true;
          }
        }
        if (local_changed) {
            clusters[nearest_partition].push_back(clusters[cluster][i]);
            clusters[cluster].erase(clusters[cluster].begin() + i);
        }
      }
    }
  }

  for (size_t i = 0; i < index->dop; ++i)
    index->partitions[i].clear();

  for (cluster = 0; cluster < num_clusters; ++cluster)
    for (size_t i = 0; i < clusters[cluster].size(); ++i)
      index->partitions[i % index->dop].push_back(clusters[cluster][i]);
}

void load_partitions(KrakenIndex* index) {
  index->va_files = new VAFile*[index->dop];
  index->vafile_points = new std::vector<std::vector<float> >[index->dop];

  // load data points
  for (size_t i = 0; i < index->dop; ++i) {
    index->vafile_points[i] = std::vector<std::vector<float> >(index->partitions[i].size(),
							       std::vector<float>(index->dop));
    for (size_t j = 0; j < index->partitions[i].size(); ++j)
      index->vafile_points[i][j] = index->partitions[i][j];
  }

  // create VA Files
  for (size_t i = 0; i < index->dop; ++i)
    index->va_files[i] = new VAFile(4, index->vafile_points[i]);
}

void load_vafiles(KrakenIndex* index, std::vector<std::vector<float> > vafile_points) {
  index->va_files = new VAFile*[index->dop];
  index->vafile_points = new std::vector<std::vector<float> >[index->dop];
  size_t partition_size = index->count / index->dop;

  // load data points
  for (size_t i = 0; i < index->dop; ++i) {
    index->vafile_points[i] = std::vector<std::vector<float> >(partition_size, std::vector<float>(index->dop));
    for (size_t j = 0; j < partition_size; ++j)
      index->vafile_points[i][j] = vafile_points[i*partition_size + j];
  }

  // create VA Files
  for (size_t i = 0; i < index->dop; ++i)
    index->va_files[i] = new VAFile(4, index->vafile_points[i]);
}

std::vector<uint32_t> intersect(std::vector<uint32_t> first, std::vector<uint32_t> second) {
  std::vector<uint32_t> intersection;
  std::set_intersection(first.begin(), first.end(), second.begin(),
                        second.end(), std::back_inserter(intersection));

  return intersection;
}

inline void scan_partition_vafile(int id, VAFile* vafile, std::vector<uint64_t> &results, std::vector<float> lower, std::vector<float> upper) {
  results = vafile->rangeQuery(lower, upper);
}

inline void scan_partition_vafile_simd(int id, VAFile* vafile, std::vector<uint64_t> &results, std::vector<float> lower, std::vector<float> upper) {
  results = vafile->rangeQuerySIMD(lower, upper);
}

std::vector<uint64_t> partitioned_range_vafile(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint64_t> results;
  std::vector<std::vector<uint64_t>> intermediate_results(index->dop, std::vector<uint64_t>());
  std::future<void> *futures = new std::future<void>[index->dop];

  for (size_t i = 0; i < index->dop; i++)
    futures[i] = pool->push(std::ref(scan_partition_vafile), index->va_files[i], std::ref(intermediate_results[i]), lower, upper);

  for (size_t i = 0; i < index->dop; i++) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }

  delete [] futures;

  return results;
}

std::vector<uint64_t> partitioned_range_vafile_simd(KrakenIndex* index, ctpl::thread_pool *pool, std::vector<float> lower, std::vector<float> upper) {
  std::vector<uint64_t> results;
  std::vector<std::vector<uint64_t>> intermediate_results(index->dop, std::vector<uint64_t>());
  std::future<void> *futures = new std::future<void>[index->dop];

  for (size_t i = 0; i < index->dop; i++)
    futures[i] = pool->push(std::ref(scan_partition_vafile_simd), index->va_files[i], std::ref(intermediate_results[i]), lower, upper);

  for (size_t i = 0; i < index->dop; i++) {
    futures[i].get();
    results.insert(std::end(results), std::begin(intermediate_results[i]), std::end(intermediate_results[i]));
  }

  delete [] futures;

  return results;
}
