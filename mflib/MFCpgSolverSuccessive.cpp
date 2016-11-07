#include "MFCpgSolverSuccessive.hpp"
#include "MFCpgEstimator.hpp"
#include "MFSolver.hpp"

namespace methylFlow {
    
    MFCpgSolverSuccessive::MFCpgSolverSuccessive(MFGraph* mfobj, const float length_mult) : MFSolver(mfobj),
    estimator(mfobj, &std::cout, length_mult),
    alpha_lambda(mfobj->get_graph()),
    beta_lambda(mfobj->get_graph())
    {
        estimator.computeNormalized();
    }
    
    MFCpgSolverSuccessive::~MFCpgSolverSuccessive()
    {
	}
	//L Prime F Calculation
	 float MFCpgSolverSuccessive::LPrimeOfF(int estimate_f){
		float sum = 0;
		
		for (MFCpgEstimator::CpgMap<float>::iterator it = estimator.normalized_map.begin();
             it != estimator.normalized_map.end(); ++it){
                    int pos = it->first;
					MFCpgEstimator::CpgEntry<float> entry = it->second;
					float u = entry.Cov;
					float m = entry.Meth;
					sum += signFunction(estimate_f) * (-(firstNumerator(pos, m)/MLofF(pos,estimate_f,u, m) + secondNumerator(pos,m)/ULofF(pos,estimate_f, u, m)));
                }
		return sum;
	}
	
	 double MFCpgSolverSuccessive::firstNumerator(int ell, int m){
		
		if(ell == m){
			return 1/ell;
		}
		return 0;
	}
	
	 int MFCpgSolverSuccessive::MLofF(int ell, double estimate_f, int u, int m){
		
		double sum = 0;
		for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
            ListDigraph::Node v = mfGraph.source(arc);
            if (ell == m){
				for(int j = 0; j <= u; j++){
					sum += estimate_f /lvu;
				}
			}
        }
		return sum;

	}
	
	 int MFCpgSolverSuccessive::ULofF(int ell, double estimate_f, int u, int m){
		
		double sum = 0;
		for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
            ListDigraph::Node v = mfGraph.source(arc);
            if (ell == u){
				for(int j = 0; j <= u; j++){
					sum += estimate_f /ell;
				}
			}
        }
		
		return sum;
	}
	
	
		
	 double MFCpgSolverSuccessive::secondNumerator(int ell, int u){
		
		if(ell == u){
			return 1/ell;
		}
		return 0;
	}
	
	
	 int MFCpgSolverSuccessive::signFunction(int estimate_f){
		if(estimate_f < 0){
			return -1;
		}else if (estimate_f > 0){
			return 1;
		}else 
			return 0;
	}
    
    float MFCpgSolverSuccessive::score(const float lambda)
    {
        float obj = lp->primal();
        for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
            obj -= lambda * lp->dual(rows[arc]);
        }
        return obj;
    }
    
    int MFCpgSolverSuccessive::add_cols()
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
        
        // now add alpha and beta for lambda nodes
        for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
            ListDigraph::Node v = mf->mfGraph.source(arc);
            
            alpha_lambda[v] = lp->addCol();
            lp->colLowerBound(alpha_lambda[v], 0.);
            lp->colUpperBound(alpha_lambda[v], 1.);
            
            beta_lambda[v] = lp->addCol();
            lp->colLowerBound(beta_lambda[v], 0.);
            lp->colUpperBound(beta_lambda[v], 1.);
        }
        return 0;
    }
    
    int MFCpgSolverSuccessive::make_deviance_objective(Lp::Expr &obj)
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
    
    int MFCpgSolverSuccessive::make_lambda_objective(const float lambda, Lp::Expr &obj)
    {
        const ListDigraph &mfGraph = mf->get_graph();
        
        // add terms to objective for lambda nodes
        for (ListDigraph::InArcIt arc(mfGraph, mf->get_sink()); arc != INVALID; ++arc) {
            ListDigraph::Node v = mfGraph.source(arc);
            obj += lambda * (beta_lambda[v] - alpha_lambda[v]);
        }
        return 0;
    }
    
    int MFCpgSolverSuccessive::add_constraints()
    {
        const ListDigraph &mfGraph = mf->get_graph();
        
        //add sink constraints
        for (ListDigraph::InArcIt arc(mfGraph, mf->get_sink()); arc != INVALID; ++arc) {
            ListDigraph::Node v = mfGraph.source(arc);
            
            rows[arc] = lp->addRow(scaled_length[arc] * beta_lambda[v] -
                                   scaled_length[arc] * alpha_lambda[v] - nu[v] <= -CONSISTENCY_FACTOR);
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
                    expr += beta_y[loc] - alpha_y[loc];
                    if (entry.methyl) {
                        expr += beta_m[loc] - alpha_m[loc];
                    }
                }
                expr += nu[u] - nu[v];
                rows[arc] = lp->addRow(expr <= -CONSISTENCY_FACTOR);
            }
        }
        return 0;
    }
    
    int MFCpgSolverSuccessive::modify_lambda_constraints(const float lambda)
    {
        const ListDigraph &mfGraph = mf->get_graph();
        
        for (ListDigraph::InArcIt arc(mf->mfGraph, mf->sink); arc != INVALID; ++arc) {
            ListDigraph::Node v = mfGraph.source(arc);
            Lp::Row row = rows[arc];
            lp->row(row, -lambda * beta_lambda[v] - (-lambda * alpha_lambda[v]) - nu[v] <= -CONSISTENCY_FACTOR);
        }
        return 0;
    }
    
    void MFCpgSolverSuccessive::print_primal()
    {
        for (MFCpgEstimator::CpgMap<float>::iterator it = estimator.normalized_map.begin();
             it != estimator.normalized_map.end(); ++it) {
            int pos = it->first;
            std::cout << pos << ": ay=" << lp->primal(alpha_y[pos]);
            std::cout << " by=" << lp->primal(beta_y[pos]);
            std::cout << " am=" << lp->primal(alpha_m[pos]);
            std::cout << " bm=" << lp->primal(beta_m[pos]);
            
            MFCpgEstimator::CpgEntry<float> entry = it->second;
            std::cout << " y=" << entry.Cov << " my=" << entry.Meth << std::endl;
        }
        print_nus();
    }
} // namespace methylFlow
