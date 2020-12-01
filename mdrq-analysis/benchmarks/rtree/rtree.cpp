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
#include "rtree.h"

Index* createRTreeIndex(uint dim) {
  Tools::PropertySet* ps = GetDefaults();
  Tools::Variant var;
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = RT_RTree;
  ps->setProperty("IndexType", var);
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = dim;
  ps->setProperty("Dimension", var);
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = RT_Memory;
  ps->setProperty("IndexStorageType", var);
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = 96;
  ps->setProperty("IndexCapacity", var);
  var.m_varType = Tools::VT_ULONG;
  var.m_val.ulVal = 96;
  ps->setProperty("LeafCapacity", var);

  Index* idx = new Index(*ps);
  delete ps;

  if (!idx->index().isIndexValid())
    throw "Failed to create valid index";

  return idx;
}

void addPointToRTree(Index* idx, std::vector<float> data_vector, uint dimensions, int64_t id) {
  double* coords = new double[dimensions];
  for (uint i = 0; i < dimensions; i++)
    coords[i] = (double) data_vector[i];
  uint8_t* pData = 0;
  uint32_t nDataLength = 0;

  SpatialIndex::IShape* shape = new SpatialIndex::Point(coords, dimensions);
  idx->index().insertData(nDataLength,pData,*shape,id+1);
  delete shape;
}

std::vector<uint32_t> rangeRTree(Index* idx, std::vector<float> lower, std::vector<float> upper, uint dim) {
  double* plow  = new double[dim];
  double* phigh = new double[dim];
  for (size_t i = 0; i < dim; ++i) {
    plow[i]  = lower[i];
    phigh[i] = upper[i];
  }

  // get nearest maxResults shapes form index
  CustomVisitor* visitor = new CustomVisitor;
  SpatialIndex::Region* r = new SpatialIndex::Region(plow, phigh, dim);
  idx->index().intersectsWithQuery(*r, *visitor);

  uint64_t nResultCount = visitor->GetResultCount();

  // get actual results
  std::vector<SpatialIndex::IData*>& results = visitor->GetResults();
  std::vector<uint32_t> resultIds;

  for (size_t i = 0; i < nResultCount; ++i)
    resultIds.push_back(results[i]->getIdentifier());

  delete [] plow;
  delete [] phigh;
  delete visitor;

  return resultIds;
}

std::vector<uint32_t> rangeRTreeSIMD(Index* idx, std::vector<float> lower, std::vector<float> upper, uint dim) {
  double* plow  = new double[dim];
  double* phigh = new double[dim];
  for (size_t i = 0; i < dim; ++i) {
    plow[i]  = lower[i];
    phigh[i] = upper[i];
  }

  // get nearest maxResults shapes form index
  CustomVisitor* visitor = new CustomVisitor;
  SpatialIndex::Region* r = new SpatialIndex::Region(plow, phigh, dim);
  idx->index().intersectsWithQuerySIMD(*r, *visitor);

  uint64_t nResultCount = visitor->GetResultCount();

  // get actual results
  std::vector<SpatialIndex::IData*>& results = visitor->GetResults();
  std::vector<uint32_t> resultIds;

  for (size_t i = 0; i < nResultCount; ++i)
    resultIds.push_back(results[i]->getIdentifier());

  delete [] plow;
  delete [] phigh;
  delete visitor;

  return resultIds;
}

CustomVisitor* rangeRTreeVisitor(Index* idx, std::vector<float> lower, std::vector<float> upper, uint dim) {
  double* plow  = new double[dim];
  double* phigh = new double[dim];
  for (size_t i = 0; i < dim; ++i) {
    plow[i]  = lower[i];
    phigh[i] = upper[i];
  }

  // get nearest maxResults shapes form index
  CustomVisitor* visitor = new CustomVisitor;
  SpatialIndex::Region* r = new SpatialIndex::Region(plow, phigh, dim);
  idx->index().intersectsWithQuery(*r, *visitor);

  uint64_t nResultCount = visitor->GetResultCount();

  // get actual results
  std::vector<SpatialIndex::IData*>& results = visitor->GetResults();
  std::vector<uint32_t> resultIds;

  for (size_t i = 0; i < nResultCount; ++i)
    resultIds.push_back(results[i]->getIdentifier());

  delete [] plow;
  delete [] phigh;

  return visitor;
}
