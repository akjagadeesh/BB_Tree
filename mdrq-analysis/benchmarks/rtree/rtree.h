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
#ifndef __Rtree_H
#define __Rtree_H

#include <spatialindex/capi/sidx_api.h>
#include <spatialindex/capi/sidx_impl.h>
#include <spatialindex/capi/sidx_config.h>

class CustomVisitor : public SpatialIndex::IVisitor
{
  private:
    std::vector<SpatialIndex::IData*> m_vector;
    uint64_t nResults;

  public:
    size_t m_indexIO;
    size_t m_leafIO;

    CustomVisitor() : nResults(0), m_indexIO(0), m_leafIO(0) {};
    ~CustomVisitor() {
      std::vector<SpatialIndex::IData*>::iterator it;
      for (it = m_vector.begin(); it != m_vector.end(); it++) {
        delete *it;
      }
    };

    uint64_t GetResultCount() const {
      return nResults;
    }

    std::vector<SpatialIndex::IData*>& GetResults()  {
      return m_vector;
    }

    void visitNode(const SpatialIndex::INode& n) {
      if (n.isLeaf())
        m_leafIO++;
      else
        m_indexIO++;
    }

    void visitorNode(const SpatialIndex::INode& ) {
    }

    void visitData(const SpatialIndex::IData& d) {

      SpatialIndex::IData* item = dynamic_cast<SpatialIndex::IData*>(const_cast<SpatialIndex::IData&>(d).clone());
      nResults += 1;
      m_vector.push_back(item);
    }

    void visitData(std::vector<const SpatialIndex::IData*>& ) {
    }
};

Index* createRTreeIndex(uint dim);
void addPointToRTree(Index* idx, std::vector<float> data_vector, uint dimensions, int64_t id);
std::vector<uint32_t> rangeRTree(Index* idx, std::vector<float> lower, std::vector<float> upper, uint dim);
std::vector<uint32_t> rangeRTreeSIMD(Index* idx, std::vector<float> lower, std::vector<float> upper, uint dim);
CustomVisitor* rangeRTreeVisitor(Index* idx, std::vector<float> lower, std::vector<float> upper, uint dim);

#endif
