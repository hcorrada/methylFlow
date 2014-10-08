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

  template<class T> void MFCpgEstimator::computeMap(CpgMap<T> &cpg_map, CoverageFunc<T> cov_function) 
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

      for (std::vector<MethylRead::CpgEntry>::iterator it = m->cpgs.begin(); it != m->cpgs.end(); ++it) {
	MethylRead::CpgEntry entry = *it;
	
        #ifndef NDEBUG
	std::cout << "Offset " << entry.offset << ":" << (entry.methyl ? 'M' : 'U') << std::endl;
	#endif

	if (cpg_map.count(rPos + entry.offset - 1) == 0) {
	  cpg_map[rPos + entry.offset - 1] = MFCpgEstimator::CpgEntry<T>(cov, entry.methyl ? cov : 0);
	} else {
	  MFCpgEstimator::CpgEntry<T> estimatorEntry = cpg_map[rPos + entry.offset - 1];
	  estimatorEntry.Cov += cov;
	  estimatorEntry.Meth += (entry.methyl ? cov : 0);
	  cpg_map[rPos + entry.offset - 1] = estimatorEntry;
	}
      }
    }

    for (typename CpgMap<T>::iterator it = cpg_map.begin(); it != cpg_map.end(); ++it) {
      it->second.Beta = (float) it->second.Meth / (float) it->second.Cov;
    }
  }
  
  void MFCpgEstimator::computeRaw()
  {
    #ifndef NDEBUG
    std::cout << "computing raw" << std::endl;
    #endif 

    MFGraph *graph = solver->mf;
    computeMap(raw_map, &MFGraph::coverage);

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
      std::cout << it->first << "\t" << entry.Cov << "\t" << entry.Meth << "\t" << entry.Beta << std::endl;
    }
  }
  
  void MFCpgEstimator::printRaw() 
  {
    printMap(raw_map);
  }
} // namespace methylFlow
