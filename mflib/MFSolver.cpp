#include <iostream>
#include <cmath>


#include <stdio.h>
#include <sys/time.h>


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
    
  int MFSolver::solve(const float lambda, const float length_mult, const float epsilon, const bool verbose, const bool verboseTime)
    {
        int res;
        struct timeval tvalBeforeSolve, tvalAfterSolve;
        
        if (verboseTime) {
            gettimeofday (&tvalBeforeSolve, NULL);
        }

        res = make_lp(length_mult);
        if (res) return res;
        if (verboseTime) {
            gettimeofday (&tvalAfterSolve, NULL);
            std::cout << "Time in miliseconds run make lp:\t" << "make lp" << "\t" << "make lp" << "\t" << ((tvalAfterSolve.tv_sec - tvalBeforeSolve.tv_sec)*1000  + tvalAfterSolve.tv_usec/1000) - tvalBeforeSolve.tv_usec/1000 << std::endl;
        }
        
        if (lambda >= 0.)
            return solve_for_lambda(lambda);
        
        float best_lambda;
        if (verbose) {
            std::cout << "[methylFlow] Searching for best lambda" << std::endl;
        }
        if (verboseTime) {
            gettimeofday (&tvalBeforeSolve, NULL);
        }
        res = search_lambda(epsilon, best_lambda, length_mult, verbose, verboseTime);
        if (verboseTime) {
            gettimeofday (&tvalAfterSolve, NULL);
            std::cout << "Time in miliseconds run search lambda:\t" << "search lambda" << "\t" << "search lambda" << "\t" << ((tvalAfterSolve.tv_sec - tvalBeforeSolve.tv_sec)*1000  + tvalAfterSolve.tv_usec/1000) - tvalBeforeSolve.tv_usec/1000 << std::endl;
        }
        if (res) return res;
        
        if (verboseTime) {
            gettimeofday (&tvalBeforeSolve, NULL);
        }
        int solution = solve_for_lambda(best_lambda);
        if (verboseTime) {
            gettimeofday (&tvalAfterSolve, NULL);
            std::cout << "Time in miliseconds run solveLP:\t" << "solve LP" << "\t" << "solve LP" << "\t" << ((tvalAfterSolve.tv_sec - tvalBeforeSolve.tv_sec)*1000  + tvalAfterSolve.tv_usec/1000) - tvalBeforeSolve.tv_usec/1000 << std::endl;
        }
        return solution;

        
    }
    
  int MFSolver::search_lambda(const float epsilon, float &best_lambda, const float length_mult, const bool verbose, const bool verboseTime)
    {
        const double powlimitMin = -6.0;
        const double powlimitMax = 6.0;
        
        float best_deviance;
        float current_deviance, zero_deviance;
        int res;
        float current_lambda, zero_lambda;
        float factor;
        struct timeval tvalBeforeSolve, tvalAfterSolve, tvalBeforeFor, tvalAfterFor;

     if (verboseTime) {
         gettimeofday (&tvalBeforeSolve, NULL);
     }
        best_lambda = 0.005;
        zero_lambda = 0.005;
        solve_for_lambda(zero_lambda);
        zero_deviance = this->score(zero_lambda);
        if (verboseTime) {
            gettimeofday (&tvalAfterSolve, NULL);
            std::cout << "Time in miliseconds run zero lambda:\t" << "zero_deviance=" << "\t" << zero_deviance << "\t" << ((tvalAfterSolve.tv_sec - tvalBeforeSolve.tv_sec)*1000000  + tvalAfterSolve.tv_usec/1) - tvalBeforeSolve.tv_usec/1 << std::endl;
            gettimeofday (&tvalBeforeFor, NULL);
        }

        for (double curpow = powlimitMax; curpow <= powlimitMin; curpow-=0.5) {
            current_lambda = pow(2., curpow);
            if (verboseTime) {
                gettimeofday (&tvalBeforeSolve, NULL);
            }
            res = solve_for_lambda(current_lambda);
            if (verboseTime) {
                gettimeofday (&tvalAfterSolve, NULL);
                std::cout << "Time in miliseconds run 24 search:\t" << "current_lambda=" << "\t" << current_lambda << "\t" << ((tvalAfterSolve.tv_sec - tvalBeforeSolve.tv_sec)*1000000  + tvalAfterSolve.tv_usec/1) - tvalBeforeSolve.tv_usec/1 << std::endl;
            }
            
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
                break;
            }
        }
        if (verbose) {
            std::cout << "[methylFlow] best lambda found " << best_lambda << " deviance=" << best_deviance << std::endl;
        }
        
        if (verboseTime) {
           // std::cout << "[methylFlow] best lambda found " << best_lambda << " deviance=" << best_deviance << std::endl;
            gettimeofday (&tvalAfterFor, NULL);
            std::cout << "Time in miliseconds run 24 search for best lambda:\t" << "for" << "\t" << "for" << "\t" << ((tvalAfterFor.tv_sec - tvalBeforeFor.tv_sec)*1000000  + tvalAfterFor.tv_usec/1) - tvalBeforeFor.tv_usec/1 << std::endl;

        }
        
        if (verboseTime) {
            gettimeofday (&tvalBeforeSolve, NULL);
        }

        best_lambda = std::max(0.0, best_lambda-0.00001);
        int solution = solve_for_lambda(best_lambda);
        if (verboseTime) {
            gettimeofday (&tvalAfterSolve, NULL);
            std::cout << "Time in miliseconds run best lambda:\t" << "best_lambada=" << "\t" << best_lambda << "\t" << ((tvalAfterSolve.tv_sec - tvalBeforeSolve.tv_sec)*1000000  + tvalAfterSolve.tv_usec/1) - tvalBeforeSolve.tv_usec/1 << std::endl;
        }

        return solution;
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
      status = make_deviance_objective(deviance_obj);
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
      
      //      lp->obj(obj);
      lp->max();
      return 0;
    }
    
  int MFSolver::solve_for_lambda(const float lambda)
  {
    int status;

    Lp::Expr obj; 
    status = make_lambda_objective(lambda, obj);
    obj += deviance_obj;
    lp->obj(obj);

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
    
    #ifndef NDEBUG
    print_primal();
    #endif

    return 0;
  }

  void MFSolver::print_nus()
  {
    for (ListDigraph::NodeIt node(mf->get_graph()); node != INVALID; ++node) {
      if (mf->fake[node]) continue;

      std::cout << mf->nodeName_map[node] << ": nu=" << lp->primal(nu[node]) << std::endl;
    }
  }
} // namespace methylFlow
