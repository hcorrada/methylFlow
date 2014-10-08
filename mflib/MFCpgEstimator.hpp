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

  protected:
    MFSolver *solver;

  private:
    std::map<int, int> raw_coverage_map;
    std::map<int, int> raw_methylation_map;
    std::map<int, float> estimated_coverage_map;
    std::map<int, float> estimated_methylation_map;
  };

} // namespace methylFlow

#endif // MFCPGESTIMATOR_H
