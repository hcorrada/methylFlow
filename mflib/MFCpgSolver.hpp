#include "MFSolver.hpp"
#include "MFCpgEstimator.hpp"

#ifndef MFCPGSOLVER_H
#define MFCPGSOLVER_H

namespace methylFlow {

  class MFCpgSolver : public MFSolver {
    friend class MFCpgEstimator;
    
  public:
    MFCpgSolver(MFGraph *mfobj, const float length_mult);
    ~MFCpgSolver();

  private:
    MFCpgEstimator estimator;
    std::map<int, Lp::Col> alpha_y;
    std::map<int, Lp::Col> beta_y;
    std::map<int, Lp::Col> alpha_m;
    std::map<int, Lp::Col> beta_m;
    ListDigraph::NodeMap<Lp::Col> alpha_lambda;
    ListDigraph::NodeMap<Lp::Col> beta_lambda;

  protected:
    float score(const float lambda);
    int add_cols();
    int make_lp_objective(Lp::Expr &obj);
    int add_constraints();
    int modify_lambda_constraints(const float lambda);
    void print_primal();
  };
} // namespace methylFlow
#endif // MFCPGSOLVER_H
