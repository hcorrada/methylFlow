#include "MFCpgEstimator.hpp"
#include "MFSolver.hpp"

using namespace lemon;

namespace methylFlow {

  MFCpgEstimator::MFCpgEstimator(MFSolver *obj) : solver(obj)
  {
  }

  MFCpgEstimator::~MFCpgEstimator()
  {
  }

  void MFCpgEstimator::computeRaw()
  {
    MFGraph *graph = solver->mf;
    const ListDigraph &mfGraph = solver->get_graph();

    for (ListDigraph::NodeIt node(mfGraph); node != INVALID; ++node) {
      MethylRead *m = graph->read(node);
      int rPos = m->start();
      int cov = graph->coverage(node);
    }
  }

  float MFCpgEstimator::getPctError()
  {
    return 0.1;
  }
} // namespace methylFlow
