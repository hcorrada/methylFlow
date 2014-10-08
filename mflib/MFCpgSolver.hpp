#include "MFSolver.hpp"

#ifndef MFCPGSOLVER_H
#define MFCPGSOLVER_H

namespace methylFlow {

  class MFCpgSolver : public MFSolver {
    friend class MFCpgEstimator;
    
    // get pct error for current solution
    float get_pcterror();
  };

} // namespace methylFlow
#endif // MFCPGSOLVER_H
