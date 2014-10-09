#include "MFGraph.hpp"

using namespace lemon;

#ifndef MFCPGESTIMATOR_H
#define MFCPGESTIMATOR_H

namespace methylFlow {
  class MFCpgSolver;

  class MFCpgEstimator {
    friend class MFCpgSolver;

  public:
    MFCpgEstimator(MFCpgSolver *obj, const float scale_mult);
    ~MFCpgEstimator();

    void computeRaw();
    void computeNormalized();
    void computeEstimated();

    float getPctError();
    void printRaw();
    void printEstimated();

  protected:
    MFCpgSolver *solver;

    template<class T>
    class CpgEntry {
      friend class MFCpgEstimator;
      friend class MFCpgSolver;
    public:
      CpgEntry(T cov, T meth) : Cov(cov), Meth(meth), Beta(0) {};
      CpgEntry() : Cov(0), Meth(0), Beta(0) {};
    protected:
      T Cov;
      T Meth;
      float Beta;
    };

    const float scale_mult;

    template<class T> using CpgMap = std::map<int, CpgEntry<T> >;
    
    CpgMap<int> raw_map;
    CpgMap<float> normalized_map;
    CpgMap<float> estimated_map;

  private:    
    template<class T> void printMap(CpgMap<T>);

    template<class T> using CoverageFunc = T (MFCpgEstimator::*)(MFGraph *mf, const ListDigraph::Node &);
    template<class T> void computeMap(CpgMap<T> &, CoverageFunc<T>);

    int coverage(MFGraph *mf, const ListDigraph::Node &);
    float expected_coverage(MFGraph *mf, const ListDigraph::Node &);
    float normalized_coverage(MFGraph *mf, const ListDigraph::Node &);
    float calculateError();
  };

} // namespace methylFlow

#endif // MFCPGESTIMATOR_H
