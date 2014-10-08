#include <iostream>
#include <cmath>

#include <lemon/lp.h>

#include "MFSolver.hpp"

using namespace lemon;

namespace methylFlow {
    MFSolver::MFSolver(MFGraph *mfobj) : mf(mfobj),
    nu(mfobj->get_graph()),
    rows(mfobj->get_graph()),
    scaled_length(mfobj->get_graph())
    {
    }
    
    MFSolver::~MFSolver()
    {
    }
    
  int MFSolver::solve(const float lambda, const float length_mult, const float epsilon, const bool verbose)
    {
        int res;
        res = make_lp(length_mult);
        if (res) return res;
        
        if (lambda >= 0.)
            return solve_for_lambda(lambda);
        
        float best_lambda;
        if (verbose) {
            std::cout << "[methylFlow] Searching for best lambda" << std::endl;
        }
        
        res = search_lambda(epsilon, best_lambda, length_mult, verbose);
        if (res) return res;
        
        return solve_for_lambda(best_lambda);
    }
    
  int MFSolver::search_lambda(const float epsilon, float &best_lambda, const float length_mult, const bool verbose)
    {
        const double powlimit = 6.0;
        
        float best_deviance;
        float current_deviance, zero_deviance;
        int res;
        float current_lambda, zero_lambda;
        float factor;
     
        best_lambda = 0;
        zero_lambda = 0;
        solve_for_lambda(zero_lambda);
        zero_deviance = this->score(zero_lambda);
        
        for (double curpow = -powlimit; curpow <= powlimit; curpow+=.5) {
            current_lambda = pow(2., curpow);
            res = solve_for_lambda(current_lambda);
            if (res) return res;
            
            current_deviance = this->score(current_lambda);
            #ifndef NDEBUG
                std::cout << "lam=2^" << curpow << " dev=" << current_deviance;
                std::cout << " , opt = " <<  lp->primal()  << std::endl;
            #endif
            

          //  factor = zero_deviance / current_deviance;
            if (current_deviance < 0.00001 || (factor = zero_deviance / current_deviance >= 1.0 - epsilon) ) {
                best_deviance = current_deviance;
                best_lambda = current_lambda;

            }
        }
        
        if (verbose) {
            std::cout << "[methylFlow] best lambda found " << best_lambda << " deviance=" << best_deviance << std::endl;
        }
	best_lambda = std::max(0.0, best_lambda-0.00001);
        return solve_for_lambda(best_lambda);
    }
    
    int MFSolver::make_lp(const float length_mult)
    {
      int status;
      const ListDigraph &mfGraph = mf->get_graph();
        
      lp = new Lp();
        
      // scale the lengths
      for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
	// divid int by int is not a float.
	scaled_length[arc] = float(mf->effective_length(arc)) / length_mult;
      }

	
#ifndef NDEBUG
      std::cout << "making LP on " << countNodes(mf->mfGraph) << " nodes" << std::endl;
#endif
      
      status = add_cols();
      if (status) return status;

      for (ListDigraph::NodeIt v(mfGraph); v != INVALID; ++v) {
	if (mf->fake[v]) continue;                        
            // add node's nu variable
	nu[v] = lp->addCol();            
      }
      Lp::Expr obj;
      status = make_lp_objective(obj);
      if (status) return status;

      // bound nu variable for source targets
      for (ListDigraph::OutArcIt arc(mfGraph, mf->get_source()); arc != INVALID; ++arc) {
            ListDigraph::Node v = mfGraph.target(arc);
            rows[arc] = lp->addRow(nu[v] <= -CONSISTENCY_FACTOR);
      }
#ifndef NDEBUG
      std::cout << "nu bounds added" << std::endl;
#endif
      
      status = add_constraints();
      if (status) return status;

#ifndef NDEBUG
      std::cout << "constraints added" << std::endl;
#endif
      
      lp->obj(obj);
      lp->max();
      return 0;
    }
    
  int MFSolver::solve_for_lambda(const float lambda)
  {
    int status;

    status = modify_lambda_constraints(lambda);
    if (status) return status;

#ifndef NDEBUG
    std::cout << "lambda constraints updated" << std::endl;
    std::cout << "running solver on " << countNodes(mf->get_graph()) << " nodes" << std::endl;
#endif

    lp->solve();
#ifndef NDEBUG
    std::cout << "obj = " << lp->primal() << std::endl;
    std::cout << "get last score = " << score(lambda) << std::endl;
    #endif

    if (lp->primalType() != Lp::OPTIMAL) {
      std::cout << "Did not find optimum" << std::endl;
      return -1;
    } 
    return 0;
  }
    
  int MFSolver::extract_flows()
  {
    // extract flows
    // TODO: check Lp is there and solved
    for (ListDigraph::ArcIt arc(mf->get_graph()); arc != INVALID; ++arc) {
#ifndef NDEBUG
      std::cout << "Extracting flow of arc: " << " ";
      std::cout << mf->nodeName_map[mf->get_graph().source(arc)];
      std::cout << " -> ";
      std::cout << mf->nodeName_map[mf->get_graph().target(arc)];
#endif
      Lp::Row row = rows[arc];
      mf->flow_map[arc] = lp->dual(row);
      
#ifndef NDEBUG
      std::cout << " flow: " << mf->flow_map[arc] << std::endl;
	    #endif
    }
    return 0;
  }
  
} // namespace methylFlow
