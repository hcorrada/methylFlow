#include "MFGraph.hpp"

using namespace lemon;

#ifndef MFSOLVER_H
#define MFSOLVER_H

#define CONSISTENCY_FACTOR 0

namespace methylFlow {

  class MFSolver {
  public:
    MFSolver(MFGraph *mfobj);
    ~MFSolver();
    
    // does everything
   virtual  int solve(const float lambda, const float length_mult, const float epsilon, const bool verbose, const bool verboseTime );

    // extract flows from LP solution
    int extract_flows();

    ListDigraph &get_graph();
    const ListDigraph &get_graph() const;

  protected:
    MFGraph *mf;
    Lp *lp;
    Lp::Expr deviance_obj;

    ListDigraph::NodeMap<Lp::Col> nu;
    ListDigraph::ArcMap<Lp::Row> rows;
    ListDigraph::ArcMap<float> scaled_length;

    // make the LP object
    // virtual LP creator
    int make_lp(const float length_mult);

    virtual int add_cols() =0;
    virtual int make_deviance_objective(Lp::Expr &obj) =0;
    virtual int make_lambda_objective(const float lambda, Lp::Expr &obj) =0;
    virtual int add_constraints() =0;
    virtual int modify_lambda_constraints(const float lambda) =0;
    virtual void print_primal() =0;
    void print_nus();

    // solve the optimization problem
    // lambda: penalty parameter
    int solve_for_lambda(const float lambda);

    // find the best lambda
    int search_lambda(const float epislon, float &best_lambda, const float scale_mult, const bool verbose, const bool verboseTime);
    
    // virtual score function
    virtual float score(const float lambda) =0;
  };

  inline ListDigraph &MFSolver::get_graph()
  {
    return mf->get_graph();
  }

  inline const ListDigraph &MFSolver::get_graph() const
  {
    return mf->get_graph();
  }

} // namespace methylFlow

#endif // MFSOLVER_H
