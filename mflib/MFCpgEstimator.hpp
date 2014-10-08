#include "MFGraph.hpp"

using namespace lemon;

#ifndef MFCPGESTIMATOR_H
#define MFCPGESTIMATOR_H

namespace methylFlow {
  class MFCpgEstimator {
  public:
    MFCpgEstimator(MFSolver *obj);
    ~MFCpgEstimator();

    void computeRaw();
    float getPctError();
    void printRaw();

  protected:
    MFSolver *solver;

  private:    
    template<class T>
    class CpgEntry {
      friend class MFCpgEstimator;
    protected:
      T Cov;
      T Meth;
      float Beta;
    };

    template<class T> using CpgMap = std::map<int, CpgEntry<T> >;
    
    CpgMap<int> raw_map;
    CpgMap<float> estimated_map;

    template<class T> void printMap(CpgMap<T>);

    template<class T> using CoverageFunc = T &(MFGraph::*)(const ListDigraph::Node &);
    template<class T> void computeMap(CpgMap<T>, CoverageFunc<T>);
  };

} // namespace methylFlow

#endif // MFCPGESTIMATOR_H
