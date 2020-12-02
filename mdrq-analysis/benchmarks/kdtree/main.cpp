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
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

#include "kraken.h"

#define FEATUREVECTORS_FILE "../../1000genomes_import/chr22_feature.vectors"
#define GENES_FILE "../../1000genomes_import/genes.txt"
#define POWER_DATASET_FILE "../../power_import/DEBS2012-ChallengeData.txt"

static double gettime(void) {
  struct timeval now_tv;
  gettimeofday (&now_tv,NULL);
  return ((double)now_tv.tv_sec) + ((double)now_tv.tv_usec) / 1000000.0;
}

//ADDED--
// determines the statistical average (mean) of a given sequence of doubles.
static double getaverage(double* runtimes, size_t n) {
  double sum = 0.0;

  for (size_t i = 0; i < n; ++i) {
    sum += runtimes[i];
  }

  return (sum / n);
}

// determines the standard deviation of a given sequence of doubles.
static double getstddev(double* runtimes, size_t n) {
  double avg = 0.0;
  double std_dev = 0.0;

  for (size_t i = 0; i < n; ++i) {
    avg += runtimes[i];
  }
  avg = (avg / n);

  for (size_t i = 0; i < n; ++i) {
    std_dev += (runtimes[i] - avg) * (runtimes[i] - avg);
  }
  std_dev = std_dev / n;
  std_dev = sqrt(std_dev);

  return std_dev;
}
//--ADDED

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " num_elements num_dimensions distribution(0=normal, 1=clustered, 2=uniform, 3=gmrqb, 4=power)" << std::endl;
    return 1;
  }

  size_t n = atoi(argv[1]);
  size_t m = atoi(argv[2]);
  int o = 1;
  float selectivity = 0.5;

  std::cout << "INFO: " << n << " vectors, " << m << " dimensions." << std::endl;

  if (atoi(argv[3]) == 2 && argc == 5) {
    selectivity = atof(argv[4]);
    std::cout << "Query Selectivity per Dimension: " << selectivity << std::endl;
  }

  std::vector<std::vector<float> > kdtree_points(n, std::vector<float>(m));
  ctpl::thread_pool *vscan_pool;
  KrakenIndex* index;
  if (argc == 6) {
	  std::cout << "THREADS: " << argv[5] << std::endl;
	  vscan_pool = new ctpl::thread_pool(atoi(argv[5]));
	  index = create_kraken(m, atoi(argv[5]));
  } else {
          vscan_pool = new ctpl::thread_pool(std::thread::hardware_concurrency());
	  index = create_kraken(m, std::thread::hardware_concurrency());
  }

  if (atoi(argv[3]) == 3) {
    size_t i = 0;
    std::ifstream feature_vectors(FEATUREVECTORS_FILE);
    std::string line;
    std::string token;

    while (std::getline(feature_vectors, line) && i < n) {
      std::vector<float> data_point(m);
      std::vector<std::string> line_tokens;
      std::istringstream iss(line);
      while(std::getline(iss, token, ' '))
        line_tokens.push_back(token);
      // parse dimensions
      for (size_t j = 0; j < m; ++j)
        data_point[j] = stof(line_tokens[j]);
      kdtree_points[i++] = data_point;
    }
    std::cout << "Loaded tuples from 1000 Genomes Project." << std::endl;
  } else if (atoi(argv[3]) == 4) {
    size_t i = 0;
    std::ifstream tuples(POWER_DATASET_FILE);
    std::string line;
    std::string token;

    while (std::getline(tuples, line) && i < n) {
      std::vector<float> data_point(3);
      std::vector<std::string> line_tokens;
      std::istringstream iss(line);
      while(std::getline(iss, token, '\t'))
        line_tokens.push_back(token);
      // parse attributes
      for (size_t j = 1; j < 4; ++j)
        data_point[j-1] = stof(line_tokens[j]);
      kdtree_points[i++] = data_point;
    }
    std::cout << "Loaded " << n << " tuples from the POWER data set." << std::endl;
  } else if(atoi(argv[3]) == 1) {
    if (argc != 5) {
      std::cout << "ERROR: Please specify a data source for the subclustered experiment!" << std::endl;
      return 0;
    }
    size_t i = 0;
    std::ifstream tuples(argv[4]);
    std::string line;
    std::string token;

    while (std::getline(tuples, line) && i < n) {
      std::vector<float> data_point(m);
      std::vector<std::string> line_tokens;
      std::istringstream iss(line);
      while(std::getline(iss, token, ','))
        line_tokens.push_back(token);
      // parse attributes
      for (size_t j = 0; j < m; ++j)
        data_point[j] = stof(line_tokens[j]);
      kdtree_points[i++] = data_point;
    }
    std::cout << "Loaded " << n << " tuples from '" << argv[4] << "'." << std::endl;
  } else {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> nd(0.5,0.5);
    std::uniform_real_distribution<double> ud(0,o);

    for (size_t i = 0; i < n; ++i) {
      std::vector<float> data_point(m);
      for (size_t j = 0; j < m; ++j) {
        if (atoi(argv[3]) == 0)
          data_point[j] = (float) nd(gen);
        else
          data_point[j] = (float) ud(gen); //ADDED
      }
      kdtree_points[i] = data_point;
    }
  }

  index->count = n;

  // random insertion order of data points
  std::random_shuffle(kdtree_points.begin(), kdtree_points.end());

  double start = gettime();
  //ADDED--
  double* runtimes;
  Tree* kd_tree = new Tree();

  std::cout << "KDTree [inserts]" << std::endl;
  runtimes = new double[n];
  for (size_t i = 0; i < n; ++i) {
    start = gettime();
    kd_tree->insertObject(kdtree_points[i], (i+1));
    runtimes[i] = (gettime() - start) * 1000000;
    insert(index, kdtree_points[i]); //used for parallelization
  }
  load_partitions(index); //used for parallelization
  std::cout << "Mean: " << getaverage(runtimes, n) << " Standard Deviation: " << getstddev(runtimes, n) << std::endl;
  delete runtimes;

  std::cout << "KDTree [point queries]" << std::endl;
  runtimes = new double[n];
  for (size_t i = 0; i < n; ++i) {
    start = gettime();
    assert((i+1) == kd_tree->searchObject(kdtree_points[i]));
    runtimes[i] = (gettime() - start) * 1000000;
  }
  std::cout << "Mean: " << getaverage(runtimes, n) << " Standard Deviation: " << getstddev(runtimes, n) << std::endl;
  delete runtimes;
  //--ADDED

  int avg_result_size = 0;
  size_t repeat = 1;
  // Generate rq queries
  int rq = 1000;
  std::vector<std::vector<float> > lb_queries(rq, std::vector<float>(m, std::numeric_limits<float>::min()));
  std::vector<std::vector<float> > ub_queries(rq, std::vector<float>(m, std::numeric_limits<float>::max()));
  if (atoi(argv[3]) == 3) {
    rq = 100;
    size_t i = 0;
    std::ifstream queries(argv[4]);
    std::cout << "Loading queries from '" << argv[4] << "'" << std::endl;
    std::string line;
    std::string token;

    while (std::getline(queries, line) && i < rq) {
      std::vector<std::string> line_tokens;
      std::istringstream iss(line);
      while(std::getline(iss, token, ' '))
        line_tokens.push_back(token);
      for (size_t j = 0; j < m; ++j) {
        if (line_tokens[j] != "min")
          lb_queries[i][j] = (float) std::stof(line_tokens[j]);
      }
      std::getline(queries, line);
      std::vector<std::string> line_tokens2;
      std::istringstream iss2(line);
      while(std::getline(iss2, token, ' '))
        line_tokens2.push_back(token);
      for (size_t j = 0; j < m; ++j) {
        if (line_tokens2[j] != "max")
          ub_queries[i][j] = (float) std::stof(line_tokens2[j]);
      }
      i++;
    }
  } else {
    for (size_t i = 0; i < rq; ++i) {
	  int first  = rand() % n;
	  int second = rand() % n;
	  for (size_t j = 0; j < m; ++j) {
	  	lb_queries[i][j] = std::min(kdtree_points[first][j], kdtree_points[second][j]);
		ub_queries[i][j] = std::max(kdtree_points[first][j], kdtree_points[second][j]);
		if (argc == 5 && atoi(argv[3]) == 2) {
		  lb_queries[i][j] = (float) ((rand() % (o * 10000)) / 1000000.0);
		  ub_queries[i][j] = lb_queries[i][j] + o*selectivity;
                }
	  }
    }
  }

  //ADDED--
  std::cout << "KDTree [range queries]" << std::endl;
  runtimes = new double[rq];
  for (size_t i = 0; i < rq; ++i) {
    start = gettime();
    std::vector<uint32_t> results =  kd_tree->rangeSearch(lb_queries[i],ub_queries[i]);
    runtimes[i] = (gettime() - start) * 1000;
    std::cout << i << std::endl;
  }
  std::cout << "Mean: " << getaverage(runtimes, rq) << " Standard Deviation: " << getstddev(runtimes, rq) << std::endl;
  delete runtimes;
  //--ADDED

  //TODO parallel ranged queries

  avg_result_size = 0;
  start = gettime();
  for (size_t r = 0; r < repeat; ++r) {
    for (size_t i = 0; i < rq; ++i) {
      avg_result_size += partitioned_range_kdtree_simd(index, vscan_pool, lb_queries[i], ub_queries[i]).size();
    }
  }
  printf("MDRQ Throughput (multi-threaded/Vertical Partitioning/SIMD): %f ops/s [avg result size: %f].\n",
		(float) ((repeat*rq) / (gettime() - start)),
		(float) (avg_result_size / (float) (repeat*rq)));

  delete index, kd_tree, vscan_pool; //TODO: deletes

  return 0;
}
