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
#include <cstring>
#include <sys/time.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

#define FEATUREVECTORS_FILE "1000genomes_import/chr22_feature.vectors"
#define GENES_FILE "1000genomes_import/genes.txt"

static double gettime(void) {
  struct timeval now_tv;
  gettimeofday (&now_tv,NULL);
  return ((double)now_tv.tv_sec) + ((double)now_tv.tv_usec) / 1000000.0;
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " num_elements num_dimensions distribution(0=normal, 1=multimodal, 2=uniform, 3=1000genomes)" << std::endl;
    return 1;
  }

  size_t  n = atoi(argv[1]);
  size_t  m = atoi(argv[2]);
  float   o = 1.0;
  float sel = 10.0;

  std::vector< std::vector<float> > data_points(n, std::vector<float>(m));

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
      data_points[i++] = data_point;
    }
  } else {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> nd(0.5,0.5);
    std::uniform_real_distribution<double> ud(0,o);
    // create multi-modal distribution
    size_t distribution_modals = 3;
    std::vector<std::normal_distribution<double>> mmd(distribution_modals);
    //  for (size_t i = 0; i < distribution_modals; ++i)
    mmd[0] = std::normal_distribution<double>(0.2, 0.2);
    mmd[1] = std::normal_distribution<double>(0.4, 0.2);
    mmd[2] = std::normal_distribution<double>(0.6, 0.2);

    for (size_t i = 0; i < n; ++i) {
      std::vector<float> data_point(m);
      for (size_t j = 0; j < m; ++j) {
        if (atoi(argv[3]) == 0)
          data_point[j] = (float) nd(gen);
        else if (atoi(argv[3]) == 1)
          data_point[j] = (float) mmd[rand() % distribution_modals](gen);
        else
          data_point[j] = (float) ((rand() % ((int) o * 1000000)) / 1000000.0);
        //row_index->dims[i].push_back(data_point[j]);
      }
      data_points[i] = data_point;
    }
  }

  // Generate rq queries
  int rq = 100;
  std::vector<std::vector<float> > lb_queries(rq, std::vector<float>(m, std::numeric_limits<float>::min()));
  std::vector<std::vector<float> > ub_queries(rq, std::vector<float>(m, std::numeric_limits<float>::max()));
  if (atoi(argv[3]) == 3) {
    size_t i = 0;
    std::ifstream genes(GENES_FILE);
    std::string line;
    std::string token;

    while (std::getline(genes, line) && i < rq) {
      std::vector<std::string> line_tokens;
      std::istringstream iss(line);
      while(std::getline(iss, token, '\t'))
        line_tokens.push_back(token);
      // Query 2 (chromosome and position)
      lb_queries[i][5] = (float) 22;
      ub_queries[i][5] = (float) 22;
      lb_queries[i][6] = (float) std::stof(line_tokens[4]) - 100000.0;
      ub_queries[i][6] = (float) std::stof(line_tokens[5]) + 200000.0;
      int query_type = 7;//rand() % 8;
      int rand_point = rand() % n;

      switch (query_type) {
        // Query 2
        case 1:
          // qual (create range around a certain qual found in the data set)
          lb_queries[i][8] = data_points[rand_point][8] * 0.5;
          ub_queries[i][8] = lb_queries[i][8] * 3;
          // depth (create range around a certain depth found in the data set)
          lb_queries[i][9] = data_points[rand_point][9] * 0.5;
          ub_queries[i][9] = lb_queries[i][9] * 3;
          // allele freq (create range using a certain allele_freq found in the data set)
          lb_queries[i][10] = data_points[rand_point][10];
          ub_queries[i][10] = lb_queries[i][10] + 0.3;
          break;
        // Query 3
        case 2:
          // gender
          lb_queries[i][2] = data_points[rand_point][2];
          ub_queries[i][2] = lb_queries[i][2];
          break;
        // Query 4 
        case 3:
          // gender
          lb_queries[i][2] = data_points[rand_point][2];
          ub_queries[i][2] = lb_queries[i][2];
          // population
          lb_queries[i][1] = data_points[rand_point][1];
          ub_queries[i][1] = lb_queries[i][1];
          break;
        // Query 5
        case 4:
          // gender
          lb_queries[i][2] = data_points[rand_point][2];
          ub_queries[i][2] = lb_queries[i][2];
          // population
          lb_queries[i][1] = data_points[rand_point][1];
          ub_queries[i][1] = lb_queries[i][1];
          // relationship
          lb_queries[i][4] = data_points[rand_point][4];
          ub_queries[i][4] = lb_queries[i][4];
          break;
        // Query 6
        case 5:
          // gender
          lb_queries[i][2] = data_points[rand_point][2];
          ub_queries[i][2] = lb_queries[i][2];
          // population
          lb_queries[i][1] = data_points[rand_point][1];
          ub_queries[i][1] = lb_queries[i][1];
          // relationship
          lb_queries[i][4] = data_points[rand_point][4];
          ub_queries[i][4] = lb_queries[i][4];
          // family_id (create range using a certain family_id found in the data set)
          lb_queries[i][3] = data_points[rand_point][3] * 0.5;
          ub_queries[i][3] = lb_queries[i][3] * 3;
          break;
        // Query 7
        case 6:
          // gender
          lb_queries[i][2] = data_points[rand_point][2];
          ub_queries[i][2] = lb_queries[i][2];
          // population
          lb_queries[i][1] = data_points[rand_point][1];
          ub_queries[i][1] = lb_queries[i][1];
          // relationship
          lb_queries[i][4] = data_points[rand_point][4];
          ub_queries[i][4] = lb_queries[i][4];
          // family_id (create range using a certain family_id found in the data set)
          lb_queries[i][3] = data_points[rand_point][3] * 0.5;
          ub_queries[i][3] = lb_queries[i][3] * 3;
          // mutation_id (create range using a certain mutation_id found in the data set)
          lb_queries[i][7] = data_points[rand_point][7] * 0.5;
          ub_queries[i][7] = lb_queries[i][7] * 3;
          break;
        // Query 8
        case 7:
          // sample_id
          lb_queries[i][0] = data_points[rand_point][0];
          ub_queries[i][0] = lb_queries[i][0];
          // population
          lb_queries[i][1] = data_points[rand_point][1];
          ub_queries[i][1] = lb_queries[i][1];
          // gender
          lb_queries[i][2] = data_points[rand_point][2];
          ub_queries[i][2] = lb_queries[i][2];
          // family_id (create range using a certain family_id found in the data set)
          lb_queries[i][3] = data_points[rand_point][3] * 0.5;
          ub_queries[i][3] = lb_queries[i][3] * 3;
          // relationship
          lb_queries[i][4] = data_points[rand_point][4];
          ub_queries[i][4] = lb_queries[i][4];
          // mutation_id (create range using a certain mutation_id found in the data set)
          lb_queries[i][7] = data_points[rand_point][7] * 0.5;
          ub_queries[i][7] = lb_queries[i][7] * 3;
          // qual (create range around a certain qual found in the data set)
          lb_queries[i][8] = data_points[rand_point][8] * 0.5;
          ub_queries[i][8] = lb_queries[i][8] * 3;
          // depth (create range around a certain depth found in the data set)
          lb_queries[i][9] = data_points[rand_point][9] * 0.5;
          ub_queries[i][9] = lb_queries[i][9] * 3;
          // allele freq (create range using a certain allele_freq found in the data set)
          lb_queries[i][10] = data_points[rand_point][10];
          ub_queries[i][10] = lb_queries[i][10] + 0.3;
          // allele_count (create range using a certain allele_count found in the data set)
          lb_queries[i][11] = data_points[rand_point][11] * 0.5;
          ub_queries[i][11] = lb_queries[i][11] * 3;
          // filter
          lb_queries[i][12] = data_points[rand_point][12];
          ub_queries[i][12] = lb_queries[i][12];
          // ref_base
          lb_queries[i][13] = data_points[rand_point][13];
          ub_queries[i][13] = lb_queries[i][13];
          // alt_base
          lb_queries[i][14] = data_points[rand_point][14];
          ub_queries[i][14] = lb_queries[i][14];
          // variant_type
          lb_queries[i][15] = data_points[rand_point][15];
          ub_queries[i][15] = lb_queries[i][15];
          // ancestral_allele
          lb_queries[i][16] = data_points[rand_point][16];
          ub_queries[i][16] = lb_queries[i][16];
          // genotypegenotype
          lb_queries[i][17] = data_points[rand_point][17];
          ub_queries[i][17] = lb_queries[i][17];
          // reference genome
          lb_queries[i][18] = data_points[rand_point][18];
          ub_queries[i][18] = lb_queries[i][18];
          break;
        // Query 1
        case 0:
        default:
          break;
      }
      i++;
    }
  } else {
    for (size_t i = 0; i < rq; ++i) {
      int el = rand() % n;
      for (size_t j = 0; j < m; ++j) {
        lb_queries[i][j] = data_points[el][j];
        ub_queries[i][j] = data_points[el][j] + o/sel;
      }
    }
  }

  // Print generated queries
  for (size_t i = 0; i < rq; ++i) {
    if (lb_queries[i][0] == std::numeric_limits<float>::min())
      std::cout << "min";
    else
      std::cout << lb_queries[i][0];
    for (size_t j = 1; j < m; ++j) {
      if (lb_queries[i][j] == std::numeric_limits<float>::min())
        std::cout << " min";
      else
        std::cout << " " << lb_queries[i][j];
    }
    std::cout << std::endl;

    if (ub_queries[i][0] == std::numeric_limits<float>::max())
      std::cout << "max";
    else
      std::cout << ub_queries[i][0];
    for (size_t j = 1; j < m; ++j) {
      if (ub_queries[i][j] == std::numeric_limits<float>::max())
        std::cout << " max";
      else
        std::cout << " " << ub_queries[i][j];
    }
    std::cout << std::endl;
  }

  return 0;
}
