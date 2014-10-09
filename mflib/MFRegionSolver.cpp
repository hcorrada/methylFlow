#include "MFRegionSolver.hpp"
#include "MFGraph.hpp"

namespace methylFlow {
  
  MFRegionSolver::MFRegionSolver(MFGraph *mfobj) : MFSolver(mfobj),
						   alpha(mfobj->get_graph()),
						   beta(mfobj->get_graph())
  {
  }

  MFRegionSolver::~MFRegionSolver()
  {
  }

  int MFRegionSolver::add_cols()
  {
    const ListDigraph &mfGraph = mf->get_graph();
    for (ListDigraph::NodeIt v(mfGraph); v != INVALID; ++v) {
      if (mf->fake[v]) continue;
      
      // add node's alpha variable and bounds
      alpha[v] = lp->addCol();
      lp->colLowerBound(alpha[v], 0.);
      lp->colUpperBound(alpha[v], 1.);
            
      // add node's beta variable and bounds
      beta[v] = lp->addCol();
      lp->colLowerBound(beta[v], 0.);
      lp->colUpperBound(beta[v], 1.);
    }
    return 0;
  }

  int MFRegionSolver::make_lp_objective(Lp::Expr &obj)
  {
    const ListDigraph &mfGraph = mf->get_graph();
    for (ListDigraph::NodeIt v(mfGraph); v != INVALID; ++v) {
      #ifndef NDEBUG
      std::cout << "Processing node " << mf->nodeName_map[v] << std::endl;
      #endif
      
      if (mf->fake[v]) continue;
      obj += mf->normalized_coverage(v) * (beta[v] - alpha[v]);
    }
    #ifndef NDEBUG
    std::cout << "obj added" << std::endl;
    #endif
    return 0;
  }

  int MFRegionSolver::add_constraints()
  {
    const ListDigraph &mfGraph = mf->get_graph();

    // add sink constraints
    for (ListDigraph::InArcIt arc(mfGraph, mf->get_sink()); arc != INVALID; ++arc) {
      ListDigraph::Node v = mfGraph.source(arc);
      rows[arc] = lp->addRow(scaled_length[arc] * beta[v] -
			     scaled_length[arc] * alpha[v] - nu[v] <= -CONSISTENCY_FACTOR);
    }
#ifndef NDEBUG
    std::cout << "sink constraints added" << std::endl;
#endif
    
    // add remaining constraints (if not childless)
    for (IterableBoolMap<ListDigraph, ListDigraph::Node>::FalseIt v(mf->fake); v != INVALID; ++v) {
      if (mf->childless[v]) continue;
      
      for (ListDigraph::OutArcIt arc(mfGraph, v); arc != INVALID; ++arc) {
	ListDigraph::Node u = mfGraph.target(arc);
	if (u == INVALID) {
	  std::cerr << "error getting target from arc" << std::endl;
	  return -1;
	}
	rows[arc] = lp->addRow(scaled_length[arc] * beta[v] -
			       scaled_length[arc] * alpha[v] - nu[v] + nu[u] <= -CONSISTENCY_FACTOR);	
      }
    }
    return 0;
  }

  int MFRegionSolver::modify_lambda_constraints(const float lambda)
  {
    const ListDigraph &mfGraph = mf->get_graph();
    // modify lambda constraints
    for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
      ListDigraph::Node v = mf->mfGraph.source(arc);
      Lp::Row row = rows[arc];
      lp->row(row, -lambda * beta[v] - (-lambda * alpha[v]) - nu[v] <= -CONSISTENCY_FACTOR);
    }
    return 0;
  }

  float MFRegionSolver::score(const float lambda)
  {
    float obj = lp->primal();
    for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
      obj -= lambda * lp->dual(rows[arc]);
    }
    return obj;
  }

  void MFRegionSolver::print_primal()
  {
    for(ListDigraph::NodeIt node(mf->mfGraph); node != INVALID; ++node) {
      if (mf->fake[node]) continue;

      std::cout << mf->nodeName_map[node] << ": a=" << lp->primal(alpha[node]);
      std::cout << " b=" << lp->primal(beta[node]);
      std::cout << " y=" << mf->normalized_coverage(node) << std::endl;
    }
    print_nus();
  }
}
