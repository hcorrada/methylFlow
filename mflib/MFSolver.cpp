#include <iostream>
#include <lemon/lp.h>

#include "MFSolver.hpp"

using namespace lemon;

namespace methylFlow {
  MFSolver::MFSolver(MFGraph *mfobj) : mf(mfobj),
				       alpha(mfobj->get_graph()),
				       beta(mfobj->get_graph()),
				       nu(mfobj->get_graph()),
				       rows(mfobj->get_graph()),
				       scaled_length(mfobj->get_graph())
  {
  }

  MFSolver::~MFSolver()
  {
  }

  int MFSolver::solve(const float lambda, const float length_mult)
  {
    int res;
    res = make_lp(length_mult);
    if (res) return res;

    if (lambda >= 0.)
      return solve_for_lambda(lambda);

    float best_lambda;
    res = search_lambda(.1, best_lambda);
    if (res) return res;

    return solve_for_lambda(best_lambda);
  }

  int MFSolver::search_lambda(const float epsilon, float &best_lambda)
  {
    best_lambda = 1.;
    return solve_for_lambda(best_lambda);
  }

  int MFSolver::make_lp(const float length_mult)
  {
    const ListDigraph &mfGraph = mf->get_graph();

    lp = new Lp();

  // scale the lengths
    for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
      scaled_length[arc] = mf->effective_length(arc) / length_mult;
    }

  Lp::Expr obj;

  #ifndef NDEBUG
  std::cout << "making LP on " << countNodes(mfGraph) << " nodes" << std::endl;
  #endif
  
  for (ListDigraph::NodeIt v(mfGraph); v != INVALID; ++v) {
    #ifndef NDEBUG
    std::cout << "Processing node " << nodeName_map[v] << std::endl;
    #endif

    if (mf->fake[v]) continue;

    // add node's alpha variable and bounds
    alpha[v] = lp->addCol();
    lp->colLowerBound(alpha[v], 0.);
    lp->colUpperBound(alpha[v], 1.);

    // add node's beta variable and bounds
    beta[v] = lp->addCol();
    lp->colLowerBound(beta[v], 0.);
    lp->colUpperBound(beta[v], 1.);

    // add node's nu variable
    nu[v] = lp->addCol();

    #ifndef NDEBUG
    std::cout << "LP vars added" << std::endl;
    #endif

    // add node's term in objective 
    obj += mf->normalized_coverage(v) * (beta[v] - alpha[v]);

    #ifndef NDEBUG
    std::cout << "obj added" << std::endl;
    #endif
  }

  // bound nu variable for source targets
  for (ListDigraph::OutArcIt arc(mfGraph, mf->get_source()); arc != INVALID; ++arc) {
    ListDigraph::Node v = mfGraph.target(arc);
    rows[arc] = lp->addRow(nu[v] <= 0);
  }
  #ifndef NDEBUG
  std::cout << "nu bounds added" << std::endl;
  #endif

  // add sink constraints
  for (ListDigraph::InArcIt arc(mfGraph, mf->get_sink()); arc != INVALID; ++arc) {
    ListDigraph::Node v = mfGraph.source(arc);
    rows[arc] = lp->addRow(scaled_length[arc] * beta[v] - 
			  scaled_length[arc] * alpha[v] - nu[v] <= 0);
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
	std::cout << "error getting target from arc" << std::endl;
	return -1;
      }
      rows[arc] = lp->addRow(scaled_length[arc] * beta[v] - 
			    scaled_length[arc] * alpha[v] - nu[v] + nu[u] <= 0);
    }
    #ifndef NDEBUG
    std::cout << "constraints added" << std::endl;
    #endif
  }

  lp->obj(obj);
  lp->max();
  return 0;
  }

  float MFSolver::get_deviance(const float lambda)
  {
    float obj = lp->primal();
    for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
      obj -= lambda * lp->dual(rows[arc]);
    }
    return obj;
  }

  int MFSolver::solve_for_lambda(const float lambda)
  {
    // modify lambda constraints
    for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
      ListDigraph::Node v = mf->mfGraph.source(arc);
      Lp::Row row = rows[arc];
      lp->row(row, -lambda * beta[v] - (-lambda * alpha[v]) - nu[v] <= 0);
    }
    #ifndef NDEBUG
    std::cout << "lambda constraints updated" << std::endl;
    std::cout << "running solver on " << countNodes(mf->get_graph()) << " nodes" << std::endl;
    #endif

    
    lp->solve();

  #ifndef NDEBUG
    std::cout << "Called solver" << std::endl;
  #endif

    if (lp->primalType() != Lp::OPTIMAL) {
      std::cout << "Did not find optimum" << std::endl;
      return -1;
    } else {
      std::cout << "lambda=" << lambda << " dev= " << get_deviance(lambda) << std::endl;
    }
    return 0;
  }

  int MFSolver::extract_flows()
  {
  // extract flows
    // TODO: check Lp is there and solved
    for (ListDigraph::ArcIt arc(mf->get_graph()); arc != INVALID; ++arc) {
    #ifndef NDEBUG
    std::cout << "Extracting flow of arc: " << std::endl;
    std::cout << nodeName_map[mfGraph.source(arc)];
    std::cout << " -> ";
    std::cout << nodeName_map[mfGraph.target(arc)] << std::endl;
    #endif
    Lp::Row row = rows[arc];
    mf->flow_map[arc] = lp->dual(row);
  }
  return 0;
}

} // namespace methylFlow
