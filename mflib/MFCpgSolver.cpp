#include "MFCpgSolver.hpp"


namespace methylFlow {

  MFCpgSolver::MFCpgSolver(MFGraph* mfobj) : MFSolver(mfobj)
  {
  }

  MFCpgSolver::~MFCpgSolver()
  {
  }

  float MFCpgSolver::score(const float lambda) 
  {
    return 0.1;
  }

  int MFCpgSolver::add_cols()
  {
    return 0;
  }

  int MFCpgSolver::make_lp_objective(Lp::Expr &obj)
  {
    return 0;
  }

  int MFCpgSolver::add_constraints()
  {
    return 0;
  }

  int MFCpgSolver::modify_lambda_constraints(const float lambda)
  {
    return 0;
  }

} // namespace methylFlow
