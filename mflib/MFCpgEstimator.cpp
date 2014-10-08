#include "MFCpgEstimator.hpp"
#include "MFSolver.hpp"

using namespace lemon;

namespace methylFlow {

  MFCpgEstimator::MFCpgEstimator(MFSolver *obj) : solver(obj)
  {
  }

  MFCpgEstimator::~MFCpgEstimator()
  {
  }

  template<class T> void MFCpgEstimator::computeMap(CpgMap<T> cpg_map, CoverageFunc<T> cov_function) 
  {
    MFGraph *graph = solver->mf;
    const ListDigraph &mfGraph = solver->get_graph();
    
    for (ListDigraph::NodeIt node(mfGraph); node != INVALID; ++node) {
      MethylRead *m = graph->read(node);
      if (!m) continue;

      int rPos = m->start();
      T cov = (graph->*cov_function)(node);
      
      #ifndef NDEBUG
      std::cout << "computing raw: cov: " << cov << " " << m->getString() << std::endl;
      #endif
    }
  }
  
  void MFCpgEstimator::computeRaw()
  {
    #ifndef NDEBUG
    std::cout << "computing raw" << std::endl;
    #endif 

    MFGraph *graph = solver->mf;
    computeMap(raw_map, &MFGraph::coverage);
    /*
    
    
      
      for (std::pair<std::vector<int>::iterator, std::vector<bool>::iterator> it(m->cpgOffset.begin(), m->methyl.begin()); 
	   it.first != m->cpgOffset.end(); ++it.first, ++it.second) {

	int offset = *(it.first);
	bool meth = *(it.second);

	#ifndef NDEBUG
	std::cout << "Offset " << offset << ":" << (meth ? 'M' : 'U') << std::endl;
	#endif

	if (raw_coverage_map.count(rPos + offset - 1) == 0) {
	  raw_coverage_map[rPos + offset - 1] = cov;
	  raw_methylation_map[rPos + offset - 1] = (meth ? cov : 0);
	} else {
	  raw_coverage_map[rPos + offset - 1] += cov;
	  raw_methylation_map[rPos + offset - 1] += (meth ? cov : 0);
	}
      }
    }
    */
    #ifndef NDEBUG
    printRaw();
    #endif
  }

  float MFCpgEstimator::getPctError()
  {
    return 0.1;
  }

  template<class T> void MFCpgEstimator::printMap(CpgMap<T> map)
  {
    std::cout << "pos\tCov\tMeth\tBeta" << std::endl;
    for (typename CpgMap<T>::iterator it = map.begin(); it != map.end(); ++it) {
      CpgEntry<T> entry = it->second;
      std::cout << it->first << "\t" << entry.Cov << "\t" << entry.Meth << entry.Beta << std::endl;
    }
  }
  
  void MFCpgEstimator::printRaw() 
  {
    printMap(raw_map);
  }
} // namespace methylFlow
