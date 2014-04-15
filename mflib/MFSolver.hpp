#include "MFGraph.hpp"

using namespace lemon;

#ifndef MFSOLVER_H
#define MFSOLVER_H

namespace methylFlow {

  class MFSolver {
  public:
    MFSolver(MFGraph *mfobj);
    ~MFSolver();
    
    // does everything
    int solve(const float lambda, const float length_mult, const bool solve);

    // extract flows from LP solution
    int extract_flows();

  protected:
    MFGraph *mf;

  private:
    Lp *lp;
    ListDigraph::NodeMap<Lp::Col> alpha;
    ListDigraph::NodeMap<Lp::Col> beta;
    ListDigraph::NodeMap<Lp::Col> nu;

    
    ListDigraph::ArcMap<Lp::Row> rows;
    ListDigraph::ArcMap<float> scaled_length;

    // make the LP object
    int make_lp(const float length_mult);

    // solve the optimization problem
    // lambda: penalty parameter
    int solve_for_lambda(const float lambda);

    // find the best lambda
    int search_lambda(const float epislon, float &best_lambda, const bool verbose);

    // get deviance for current solution
    float get_deviance(const float lambda);
  };

} // namespace methylFlow

#endif // MFSOLVER_H
