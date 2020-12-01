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
#include "vafile.h"

std::vector<uint8_t> VAFile::quantize(std::vector<float> coordinates) {
  std::vector<uint8_t> approximation(dimensions);
  for (uint32_t d = 0; d < splitPartitions.size(); d++) {
    float val = coordinates[d];
    uint8_t lastBorderIndex = (uint8_t) (splitPartitions[d].size() - 1);

    // check, if value is below data grid
    if (val < splitPartitions[d][0]) {
      approximation[d] = 0;
    }
    // check, if value is above data grid
    else if (val > splitPartitions[d][lastBorderIndex]) {
      approximation[d] = lastBorderIndex;
    }
    // normal case :  search grid position
    else {
      uint8_t pos = 0;
      for (int i = 0; i < splitPartitions[d].size(); i++) {
        if (val > splitPartitions[d][i]) {
          pos++;
        } else {
          break;
        }
      }
#ifdef OUTPUT
      std::cout << pos << std::endl;
#endif

      approximation[d] = pos;
    }
  }
  return approximation;
}


int VAFile::isPowerOfTwo(unsigned int x) {
  return ((x != 0) && !(x & (x - 1)));
}

VAFile::VAFile(uint8_t power, std::vector<std::vector<float> > coordinates) {
  // check, if p is a power of 2!
  if (!isPowerOfTwo(power)) {
    std::cout << "p should be a power of 2!" << std::endl;
  }

  partitions = power;
  points = coordinates;
  dimensions = (uint32_t) coordinates[0].size();

  std::vector<std::vector<float>> part(dimensions, std::vector<float>(partitions + 1));
  splitPartitions = part;

  for (uint32_t d = 0; d < dimensions; d++) {
    std::vector<float> temp;
    for (float p : points[d]) {
      temp.push_back(p);
    }
    // sort data
    std::sort(temp.begin(), temp.end());

    // find optimal split points
    for (int b = 0; b < partitions; b++) {
      int start = (int) (b * splitPartitions.size() / (float) partitions);
      splitPartitions[d][b] = temp[start];
    }
    // add a small epsilon to make sure that the last object is included
    splitPartitions[d][partitions] = temp[temp.size() - 1] + 0.000001f;
  }

#ifdef OUTPUT
  for (int d = 0; d < splitPartitions.size(); d++) {
    for (int j = 0; j < splitPartitions[d].size(); j++) {
      std::cout << splitPartitions[d][j] << " ";
    }
    std::cout << std::endl;
  }
#endif


  // add all points to index
  for (auto point : points) {
    vaPoints.push_back(quantize(point));
  }
}

std::vector<std::vector<float>> VAFile::initializeQueryLookupTable(std::vector<float> query) {
  std::vector<std::vector<float>> lookup(
      splitPartitions.size(),
      std::vector<float>(splitPartitions[0].size()));

  for(int d = 0; d < splitPartitions.size(); d++) {
    for(int i = 0; i < splitPartitions[0].size(); i++) {
      float temp = splitPartitions[d][i] - query[d];
      lookup[d][i] = temp*temp;
    }
  }
  return lookup;
}

float VAFile::getMinDist(std::vector<uint8_t> queryApprox,
    std::vector<uint8_t> pointApprox,
    std::vector<std::vector<float>> lookup,
    float max_dist) {

  float minDist = 0;
  for (uint32_t d = 0; d < splitPartitions.size(); d++) {

    uint8_t vp = pointApprox[d];
    uint8_t qp = queryApprox[d];

    if (vp < qp) {
      minDist += lookup[d][vp];
    } else if (vp > qp) {
      minDist += lookup[d][vp - 1];
    }

    //        minDist += getPartialMinDist2(queryApprox[d], pointApprox[d], lookup[d]);
    //        if (minDist > max_dist) {
    //            break;
    //        }
  }
  return minDist;
}

float VAFile::getDist(std::vector<float> query, std::vector<float> point) {
  float dist = 0;
  for (uint32_t d = 0; d < query.size(); d++) {
    float temp = query[d] - point[d];
    dist += temp * temp;
  }
  return dist;
}


bool VAFile::equal(std::vector<float> object1, std::vector<float> object2) {
  // Check if every dimension is the same
  for (uint64_t i = 0; i < object1.size(); ++i) {
    if (object1[i] != object2[i]) {
      return false;
    }
  }
  return true;
}

bool VAFile::equal(std::vector<uint8_t> object1, std::vector<uint8_t> object2) {
  // Check if every dimension is the same
  for (uint64_t i = 0; i < object1.size(); ++i) {
    if (object1[i] != object2[i]) {
      return false;
    }
  }
  return true;
}

std::vector<uint64_t> VAFile::rangeQuery(std::vector<float> lowerBound, std::vector<float> upperBound) {
  // Filter and search paradigm, so we need a queue
  std::vector<uint64_t> result;

  // Quantize the query point to get the grid
  std::vector<uint8_t> lowerBoundApprox = quantize(lowerBound);
  std::vector<uint8_t> upperBoundApprox = quantize(upperBound);

  for (uint64_t i = 0; i < vaPoints.size(); i++) {
    bool match = true;
    for (int j = 0; j < lowerBoundApprox.size(); j++) {
      if (lowerBoundApprox[j] > vaPoints[i][j] || upperBoundApprox[j] < vaPoints[i][j]) {
        match = false; break;
      }
    }
    if (match) {
      for (int j = 0; j < lowerBound.size(); ++j) {
	if (lowerBound[j] > points[i][j] || upperBound[j] < points[i][j]) {
	  match = false; break;
        }
      }
    }
    if (match) {
      result.push_back(i);
    }
  }

  return result;
}

std::vector<uint64_t> VAFile::rangeQuerySIMD(std::vector<float> lowerBound, std::vector<float> upperBound) {
  // Filter and search paradigm, so we need a queue
  std::vector<uint64_t> result;

  // Quantize the query point to get the grid
  std::vector<uint8_t> lowerBoundApprox = quantize(lowerBound);
  std::vector<uint8_t> upperBoundApprox = quantize(upperBound);

  for (uint64_t i = 0; i < vaPoints.size(); i++) {
    bool match = true;
    // search approximations
    for (int j = 0; j < lowerBoundApprox.size(); j++) {
      if (lowerBoundApprox[j] > vaPoints[i][j] || upperBoundApprox[j] < vaPoints[i][j]) {
        match = false; break;
      }
    }
    // search actual points
    if (match) {
	    size_t j = 0;
	    __m256 lower_reg, upper_reg, search_reg, lower_res, upper_res;
	    size_t mask_lower, mask_upper, mask;
	    size_t simd_compares = ((size_t) (lowerBound.size() / 8)) * 8;
	    for (; j < simd_compares; j += 8) {
		    lower_reg = _mm256_loadu_ps(&lowerBound[j]);
		    upper_reg = _mm256_loadu_ps(&upperBound[j]);
		    search_reg = _mm256_loadu_ps(&points[i][j]);
		    lower_res  = _mm256_cmp_ps(lower_reg, search_reg, _CMP_LE_OQ);
		    mask_lower = _mm256_movemask_ps(lower_res);
		    upper_res  = _mm256_cmp_ps(upper_reg, search_reg, _CMP_GE_OQ);
		    mask_upper = _mm256_movemask_ps(upper_res);
		    mask       = mask_lower & mask_upper;
		    // mask is set to 0xFF if every SIMD lane matchs
		    if (mask < 0xFF) {
			    match = false;
			    j = lowerBound.size();
			    break;
		    }
	    }
	    for (; j < lowerBound.size(); ++j) {
		    if (lowerBound[j] > points[i][j] || upperBound[j] < points[i][j]) {
			    match = false;
			    break;
		    }
	    }
	    if (match)
		    result.push_back(i);
    }
  }

  return result;
}
