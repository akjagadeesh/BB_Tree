#ifndef VAFILE_H
#define VAFILE_H

#include <vector>
#include <bitset>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <queue>
#include <iterator>
#include <algorithm>    // std::sort
#include <cmath>
// AVX Intrinsics (SIMD)
#include <immintrin.h>

class VAFile {
  private:
    uint32_t partitions = 0;
    uint32_t dimensions = 0;
    std::vector<std::vector<float> > points;
    std::vector<std::vector<uint8_t> > vaPoints;
    std::vector<std::vector<float>> splitPartitions;

    std::vector<uint8_t> quantize(std::vector<float> coordinates);
    int isPowerOfTwo(unsigned int x);
    float getMinDist(std::vector<uint8_t> queryApprox,
        std::vector<uint8_t> pointApprox,
        std::vector<std::vector<float>> lookup,
        float max_dist);
    bool equal(std::vector<float> object1, std::vector<float> object2);
    bool equal(std::vector<uint8_t> object1, std::vector<uint8_t> object2);
    std::vector<std::vector<float>> initializeQueryLookupTable(std::vector<float> query);

  public:
    VAFile(uint8_t p, std::vector<std::vector<float> > coordinates);
    float getDist(std::vector<float> point1, std::vector<float> point2);
    std::vector<uint64_t> rangeQuery(std::vector<float> lowerBound, std::vector<float> upperBound);
    std::vector<uint64_t> rangeQuerySIMD(std::vector<float> lowerBound, std::vector<float> upperBound);
	std::vector<float> exactSearch(int i);
};
#endif
