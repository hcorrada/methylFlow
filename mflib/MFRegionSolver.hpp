#include "MFSolver.hpp"

#ifndef MFREGIONSOLVER_H
#define MFREIGONSOLVER_H

namespace methylFlow {
  class MFRegionSolver : public MFSolver {

  public:
    MFRegionSolver(MFGraph *mfobj);
    ~MFRegionSolver();

  private:
    ListDigraph::NodeMap<Lp::Col> alpha;
    ListDigraph::NodeMap<Lp::Col> beta;

  protected:
    float score(const float lambda);
    int add_cols();
    int make_deviance_objective(Lp::Expr &obj);
    int make_lambda_objective(const float lambda, Lp::Expr &obj);
    int add_constraints();
    int modify_lambda_constraints(const float lambda);
    void print_primal();
  };
} // namespace methylFlow
#endif // MFREGIONSOLVER_H
