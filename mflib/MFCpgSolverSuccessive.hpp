#include "MFSolver.hpp"
#include "MFCpgEstimator.hpp"

#ifndef MFCPGSOLVERSUCCESSIVE_H
#define MFCPGSOLVERSUCCESSIVE_H

namespace methylFlow {

  class MFCpgSolverSuccessive : public MFSolver {
    friend class MFCpgEstimator;
    
  public:
    MFCpgSolverSuccessive(MFGraph *mfobj, const float length_mult);
    ~MFCpgSolverSuccessive();

  private:
    MFCpgEstimator estimator;
    std::map<int, Lp::Col> alpha_y;
    std::map<int, Lp::Col> beta_y;
    std::map<int, Lp::Col> alpha_m;
    std::map<int, Lp::Col> beta_m;
    ListDigraph::NodeMap<Lp::Col> alpha_lambda;
    ListDigraph::NodeMap<Lp::Col> beta_lambda;
	float LPrimeOfF(int estimate_f);
	double firstNumerator(int ell, int m);
	int MLofF(int ell, double estimate_f, int u, int m);
	int ULofF(int ell, double estimate_f, int u, int m);
	double secondNumerator(int ell, int u);
	int signFunction(int estimate_f);

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
#endif // MFCPGSOLVER_H
