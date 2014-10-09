#include "MFCpgSolver.hpp"
#include "MFCpgEstimator.hpp"
#include "MFSolver.hpp"

namespace methylFlow {

  MFCpgSolver::MFCpgSolver(MFGraph* mfobj, const float length_mult) : MFSolver(mfobj), estimator(this, length_mult)
  {
    estimator.computeRaw();
    estimator.computeNormalized();
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
    for (MFCpgEstimator::CpgMap<float>::iterator it = estimator.normalized_map.begin(); 
	 it != estimator.normalized_map.end(); ++it) {
      int pos = it->first;
      alpha_y[pos] = lp->addCol();
      lp->colLowerBound(alpha_y[pos], 0.);
      lp->colUpperBound(alpha_y[pos], 1.);

      beta_y[pos] = lp->addCol();
      lp->colLowerBound(beta_y[pos], 0.);
      lp->colUpperBound(beta_y[pos], 1.);

      alpha_m[pos] = lp->addCol();
      lp->colLowerBound(alpha_m[pos], 0.);
      lp->colUpperBound(alpha_m[pos], 1.);

      beta_m[pos] = lp->addCol();
      lp->colLowerBound(beta_m[pos], 0.);
      lp->colUpperBound(beta_m[pos], 1.);
    }
    return 0;
  }

  int MFCpgSolver::make_lp_objective(Lp::Expr &obj)
  {
    for (MFCpgEstimator::CpgMap<float>::iterator it = estimator.normalized_map.begin();
	 it != estimator.normalized_map.end(); ++it) {
      int pos = it->first;
      MFCpgEstimator::CpgEntry<float> entry = it->second;

      obj += entry.Cov * ( beta_y[pos] - alpha_y[pos]);
      obj += entry.Meth * ( beta_m[pos] - alpha_m[pos]);
    }
    return 0;
  }

  int MFCpgSolver::add_constraints()
  {
    const ListDigraph &mfGraph = mf->get_graph();

    //add sink constraints
    for (ListDigraph::InArcIt arc(mfGraph, mf->get_sink()); arc != INVALID; ++arc) {
      ListDigraph::Node v = mfGraph.source(arc);

      MethylRead *m = mf->read(v);
      if (!m) continue;

      int rPos = m->start();
      Lp::Expr expr;

      for (std::vector<MethylRead::CpgEntry>::iterator it = m->cpgs.begin();
	   it != m->cpgs.end(); ++it) {
	MethylRead::CpgEntry entry = *it;
	int loc = rPos + entry.offset - 1;
	expr += scaled_length[arc] * beta_y[loc] - scaled_length[arc] * alpha_y[loc];
	if (entry.methyl) {
	  expr += scaled_length[arc] * beta_m[loc] - scaled_length[arc] * alpha_m[loc];
	}
      }
      expr -= nu[v];
      rows[arc] = lp->addRow(expr <= -CONSISTENCY_FACTOR);
    }

    //add remaining constraints (if not childless)
    for (IterableBoolMap<ListDigraph, ListDigraph::Node>::FalseIt v(mf->fake);
	 v != INVALID; ++v) {
      if (mf->childless[v]) continue;
      
      for (ListDigraph::OutArcIt arc(mfGraph, v); arc != INVALID; ++arc) {
	ListDigraph::Node u = mfGraph.target(arc);
	if (u == INVALID) {
	  std::cerr << "error getting target from arc" << std::endl;
	  return -1;
	}

	MethylRead *m = mf->read(v);
	if (!m) continue;

	int rPos = m->start();
	Lp::Expr expr;

	for (std::vector<MethylRead::CpgEntry>::iterator it = m->cpgs.begin();
	     it != m->cpgs.end(); ++it) {
	  MethylRead::CpgEntry entry = *it;
	  int loc = rPos + entry.offset - 1;
	  expr += scaled_length[arc] * beta_y[loc] - scaled_length[arc] * alpha_y[loc];
	  if (entry.methyl) {
	    expr += scaled_length[arc] * beta_m[loc] - scaled_length[arc] * alpha_m[loc];
	  }
	}
	expr += nu[u] - nu[v];
	rows[arc] = lp->addRow(expr <= -CONSISTENCY_FACTOR);
      }
    }
    return 0;
  }

  int MFCpgSolver::modify_lambda_constraints(const float lambda)
  {
    const ListDigraph &mfGraph = mf->get_graph();
    
    for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
      ListDigraph::Node v = mfGraph.source(arc);

      MethylRead *m = mf->read(v);
      if (!m) continue;

      int rPos = m->start();
      Lp::Expr expr;

      for (std::vector<MethylRead::CpgEntry>::iterator it = m->cpgs.begin();
	   it != m->cpgs.end(); ++it) {
	MethylRead::CpgEntry entry = *it;
	int loc = rPos + entry.offset - 1;
	expr += -lambda * beta_y[loc] - (-lambda) * alpha_y[loc];
	if (entry.methyl) {
	  expr += -lambda * beta_m[loc] - (-lambda) * alpha_m[loc];
	}
      }
      expr -= nu[v];
      Lp::Row row = rows[arc];
      lp->row(row, expr <= -CONSISTENCY_FACTOR);
    }
    return 0;
  }

} // namespace methylFlow
