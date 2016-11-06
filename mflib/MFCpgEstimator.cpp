#include "MFCpgEstimator.hpp"
#include "MFCpgSolver.hpp"
#include "MFCpgSolverSuccessive.hpp"

#include <cmath>

using namespace lemon;

namespace methylFlow {
    
    MFCpgEstimator::MFCpgEstimator(MFGraph *obj, std::ostream * ostream, const float scale) : graph(obj), outstream(ostream), scale_mult(scale)
    {
    }
    
    MFCpgEstimator::~MFCpgEstimator()
    {
    }
    
    std::ostream & MFCpgEstimator::getstream()
    {
        return *outstream;
    }
    
    template<class T> void MFCpgEstimator::computeMap(CpgMap<T> &cpg_map, CoverageFunc<T> cov_function)
    {
        const ListDigraph &mfGraph = graph->get_graph();
        
        for (ListDigraph::NodeIt node(mfGraph); node != INVALID; ++node) {
            MethylRead *m = graph->read(node);
            if (!m) continue;
            
            int rPos = m->start();
            T cov = (this->*cov_function)(node);
            
#ifndef NDEBUG
            std::cout << "computing cpg: cov: " << (float) cov << " " << m->getString() << std::endl;
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
    
    int MFCpgEstimator::coverage(const ListDigraph::Node &node)
    {
        return graph->coverage(node);
    }
    
    float MFCpgEstimator::expected_coverage(const ListDigraph::Node &node)
    {
        MethylRead *m = graph->read(node);
        float len = (float) m->length();
        return graph->expected_coverage(node, scale_mult) / len;
    }
    
    float MFCpgEstimator::normalized_coverage(const ListDigraph::Node &node)
    {
        MethylRead *m = graph->read(node);
        float len = (float) m->length();
        return graph->normalized_coverage(node) / len * scale_mult;
    }
    
    void MFCpgEstimator::computeRaw()
    {
#ifndef NDEBUG
        std::cout << "computing raw" << std::endl;
#endif
        
        computeMap(raw_map, &MFCpgEstimator::coverage);
        
#ifndef NDEBUG
//        printRaw("", true);
#endif
    }
    
    void MFCpgEstimator::computeEstimated()
    {
#ifndef NDEBUG
        std::cout << "computing estimated" << std::endl;
#endif
        
        estimated_map.clear();
        computeMap(estimated_map, &MFCpgEstimator::expected_coverage);
        
#ifndef NDEBUG
        printEstimated("", true);
#endif
    }
    
    void MFCpgEstimator::computeNormalized()
    {
        computeMap(normalized_map, &MFCpgEstimator::normalized_coverage);
    }
    
    float MFCpgEstimator::calculateError()
    {
        float res = 0.0;
        for (CpgMap<float>::iterator it = estimated_map.begin(); it != estimated_map.end(); ++it) {
            int loc = it->first;
            CpgEntry<float> estimatedEntry = it->second;
            
            if (raw_map.count(loc) == 0) continue;
            CpgEntry<int> rawEntry = raw_map[loc];
            res += std::abs(rawEntry.Beta - estimatedEntry.Beta);
        }
        return res;
    }
    
//    float MFCpgEstimator::getPctError()
//    {
//        computeEstimated();
//        return calculateError();
//    }
//    
    template<class T> void MFCpgEstimator::printMap(CpgMap<T> map, std::string chr, bool with_beta)
    {
        for (typename CpgMap<T>::iterator it = map.begin(); it != map.end(); ++it) {
            CpgEntry<T> entry = it->second;
            getstream() << chr << "\t" << it->first << "\t" << entry.Cov << "\t" << entry.Meth;
            if (with_beta) {
                getstream() << "\t" << entry.Beta;
            }
            getstream() << std::endl;
        }
    }
    
    void MFCpgEstimator::printRaw(std::string chr, bool with_beta)
    {
        printMap(raw_map, chr, with_beta);
    }
    
    void MFCpgEstimator::printEstimated(std::string chr, bool with_beta)
    {
        printMap(estimated_map, chr, with_beta);
    }
} // namespace methylFlow
